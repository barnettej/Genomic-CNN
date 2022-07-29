
import tensorflow as tf
import math
from tensorflow import keras
import keras_tuner
import pandas as pd
from tensorflow.keras import Input, optimizers, regularizers
from tensorflow.keras.layers import Dense, concatenate, Lambda, Conv2D, Conv1D, MaxPool1D, Dropout, Reshape, Flatten, \
    MaxPool2D, BatchNormalization, GlobalAveragePooling2D
from tensorflow.keras import Model
from tensorflow.keras import layers
from keras_tuner import BayesianOptimization, RandomSearch, Hyperband


batch = 64
epoch = 50
train_size = 6522
valid_size = 1402
test_size = 1466
trials = 200
repeats = 2

initializer = tf.keras.initializers.HeUniform()
annotations = pd.read_csv("t1d.annotations.csv")
annotations_tensor = tf.convert_to_tensor(annotations)
annotations_tensor = tf.cast(annotations_tensor, tf.float32)

train_path = "train.input.csv"
valid_path = "valid.input.csv"
test_path = "test.input.csv"

train_ds = tf.data.experimental.make_csv_dataset(
    train_path,
    batch_size=batch,
    label_name="t1d",
    shuffle=True,
    num_parallel_reads=100,
    shuffle_buffer_size = 1000,
    ignore_errors=True)

valid_ds = tf.data.experimental.make_csv_dataset(
    valid_path,
    batch_size=batch,
    label_name="t1d",
    shuffle=True,
    num_parallel_reads=100,
    shuffle_buffer_size = 1000,
    ignore_errors=True)

test_ds = tf.data.experimental.make_csv_dataset(
    test_path,
    batch_size=batch,
    label_name="t1d",
    shuffle=True,
    num_parallel_reads=100,
    shuffle_buffer_size = 1000,
    ignore_errors=True)

for feature_batch, label_batch in test_ds.take(1):

    features = []
    for key, value in feature_batch.items():
        features.append(key)


def combine_batch_samples(samples, targets):
    inp1 = []
    for k in features:
        inp1.append(samples[k])
    inp1 = tf.stack(inp1, axis=-1)

    return {'input_1': inp1}, targets


train_dataset = train_ds.map(combine_batch_samples)
train_dataset = train_dataset.take(train_size)
valid_dataset = valid_ds.map(combine_batch_samples)
valid_dataset = valid_dataset.take(valid_size)
test_dataset = test_ds.map(combine_batch_samples)
test_dataset = test_dataset.take(test_size)


class ConvBlock(tf.keras.Model):
    def __init__(self, num_filters, kernel_size, pool_width):
        super(ConvBlock, self).__init__()
        self.conv = Conv1D(num_filters, kernel_size, activation="relu", kernel_initializer=initializer, padding='same')
        self.batchnorm = BatchNormalization()
        self.maxpool = MaxPool2D(1, pool_width)

    def call(self, x):
        x = self.conv(x)
        x = self.batchnorm(x)
        x = self.maxpool(x)

        return x


class TopConvBlock(tf.keras.Model):
    def __init__(self, num_filters, kernel_size, keep_percent):
        super(TopConvBlock, self).__init__()
        self.conv = Conv1D(num_filters, kernel_size, activation="relu", kernel_initializer=initializer)
        self.batchnorm = BatchNormalization()
        self.transpose = Lambda(lambda x: tf.transpose(x, perm=[0, 1, 3, 2]))
        self.keep_percent = keep_percent

    def build(self, input_shape):
        self.variants = input_shape[2]
        print(input_shape)
        keep = math.ceil(self.variants * self.keep_percent)
        self.topk = Lambda(lambda x: tf.sort(tf.math.top_k(x, keep).indices))

    def call(self, x):
        x = self.conv(x)
        x = self.batchnorm(x)
        x = self.transpose(x)
        y = self.topk(x)

        return x, y


def build_model(hp):
    conv_blocks = 1
    dense_layers = hp.Int("dense_layers", 1, 3)

    # annotation attachment
    input = Input(shape=(len(features)), batch_size=len(label_batch))
    input_reshape = tf.expand_dims(input, axis=1)
    annotations_tile = tf.expand_dims(annotations_tensor, axis=0)
    annotations_tile = tf.tile(annotations_tile, [input.shape[0], 1, 1])
    concat = concatenate([input_reshape, annotations_tile], axis=1)
    reshape = tf.expand_dims(concat, axis=3)
    x = tf.transpose(reshape, perm=[0, 3, 2, 1])

    for k in range(conv_blocks):
        k = str(k)

        filters = hp.Int("conv_block_filters" + k, 100, 600)
        width = 1

        x = ConvBlock(filters, width, pool_width = 1)(x)

    print(x.shape)
    x = GlobalAveragePooling2D()(x)
    print(x.shape)

    for k in range(dense_layers):
        k = str(k)

        units = hp.Int("dense_" + k, 5, 400)
        L2 = hp.Float("L2_" + k, min_value = 1e-6, max_value = 0.1, sampling = "log")
        x = layers.Dense(units=units, activation="relu", kernel_regularizer = regularizers.L2(L2))(x)
        rate = hp.Float("dense_" + k + "dropout", min_value=0.0, max_value=0.8)
        x = layers.Dropout(rate)(x)

    pheno_out = Dense(1, activation="sigmoid", name="pheno_out", kernel_initializer=tf.keras.initializers.GlorotUniform())(x)

    losses = {
        "pheno_out": "binary_crossentropy",
    }

    lossWeights = {
        "pheno_out": 1.0}

    model = Model(inputs=input, outputs=[pheno_out], )

    lr = hp.Float("lr", min_value=1e-4, max_value=1e-2, sampling="log")

    model.compile(optimizer=keras.optimizers.Adam(learning_rate=lr),
                  loss=losses, loss_weights=lossWeights, metrics=["AUC", "accuracy"])

    return (model)


tuner = BayesianOptimization(
    build_model,
    objective=keras_tuner.Objective("val_auc", direction="max"),
    max_trials=trials,
    executions_per_trial=repeats,
    project_name="t1d.keras.tuner",
    overwrite=False,
)
early_stop = tf.keras.callbacks.EarlyStopping(monitor='val_loss', mode="min", patience=5)


tuner.search(train_dataset,
             validation_data=valid_dataset,
             validation_batch_size=batch,
             validation_steps=(valid_size // batch),
             epochs=epoch,
             batch_size=batch,
             steps_per_epoch=(train_size // batch),
             callbacks=[early_stop]
             )

tuner.results_summary()

