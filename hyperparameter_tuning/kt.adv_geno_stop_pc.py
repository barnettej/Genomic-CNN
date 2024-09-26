import os
import tensorflow as tf
import math
from tensorflow import keras
import tensorflow_probability as tfp
import numpy as np
from tensorflow import data
import keras_tuner
import pandas as pd
from tensorflow.keras import Input, optimizers, regularizers
from tensorflow.keras.layers import Dense, concatenate, Lambda, Multiply, Conv2D, Conv1D, Add, MaxPool1D, Dropout, Reshape, Flatten, \
    MaxPool2D, BatchNormalization, GlobalAveragePooling2D, GlobalAveragePooling1D, AveragePooling1D
from tensorflow.keras import Model
from tensorflow.keras import layers
import keras.backend as K
import datetime
from tensorflow.keras.utils import plot_model
import itertools
from keras_tuner import BayesianOptimization, RandomSearch, Hyperband
from sklearn.preprocessing import LabelEncoder
import tensorflow_addons as tfa

batch = 
epoch = 
train_size = 
valid_size = 
test_size = 
trials = 
repeats = 

initializer = tf.keras.initializers.HeUniform()
# can use usecols to call in only a subset of the snps
annotations = pd.read_csv("/home/user/user/ukbiobank/t2d/data/t2d.annocut.p0.01.train.input.cid.annotations.csv")
annotations_tensor = tf.convert_to_tensor(annotations)
annotations_tensor = tf.cast(annotations_tensor, tf.float32)

train_path = "/home/user/user/ukbiobank/t2d/data/t2d.annocut.p0.01.train.input.with_pca.scale.csv"
valid_path = "/home/user/user/ukbiobank/t2d/data/t2d.annocut.p0.01.valid.input.with_pca.scale.csv"

train_ds = pd.read_csv(train_path)
print(train_ds.head())
x_train, pheno_train, pca_train = train_ds.values[:55168,11:], train_ds.values[:55168,10], train_ds.values[:55168,:10]

x_train = tf.cast(tf.convert_to_tensor(x_train), tf.float32)
pca_train = tf.cast(tf.convert_to_tensor(pca_train), tf.float32)
pheno_train = LabelEncoder().fit_transform(pheno_train)

valid_ds = pd.read_csv(valid_path)
x_valid, pheno_valid, pca_valid = valid_ds.values[:11712,11:], valid_ds.values[:11712,10], valid_ds.values[:11712,:10]
x_valid = tf.cast(tf.convert_to_tensor(x_valid), tf.float32)
pca_valid = tf.cast(tf.convert_to_tensor(pca_valid), tf.float32)
pheno_valid = LabelEncoder().fit_transform(pheno_valid)



print("x shape:",x_train.shape)
print("pheno shape:", pheno_train.shape)
print("pca shape:", pca_train.shape)


class ConvBlock(tf.keras.Model):
    def __init__(self, num_filters, kernel_size, pool_width):
        super(ConvBlock, self).__init__()
        self.conv = Conv1D(num_filters, kernel_size, activation="relu", kernel_initializer=initializer, padding='valid')
        self.batchnorm = BatchNormalization()
        # self.maxpool = MaxPool1D(pool_width)
        self.avgpool = AveragePooling1D(pool_width)

    def call(self, x):
        x = self.conv(x)
        x = self.batchnorm(x)
        # x = self.maxpool(x)
        x = self.avgpool(x)

        return x

@tf.custom_gradient
def gradient_reverse(x, lamda=1.0):
    y = tf.identity(x)

    def grad(dy):
        return lamda * -dy, None

    return y, grad

class GradientReversalLayer(tf.keras.layers.Layer):
    def __init__(self):
        super().__init__()

    def call(selfself, x, lamda=1.0):
        return gradient_reverse(x, lamda)

def build_model(hp):


    # input genotype data
    input = Input(shape=(x_train.shape[1]), batch_size=batch)
    geno = tf.expand_dims(input, axis=1)


    y = np.ones([input.shape[0], 1, 1]) * annotations_tensor
    # switch dimensions of the last 2 dimensions, because conv1d layers create a kernel that uses the last dimension and then runs over the second to last dimension
    y = tf.transpose(y, perm=[0, 2, 1])

    geno = tf.transpose(geno, perm = [0,2,1])
    x = np.ones([1, 1, y.shape[2]]) * geno
    x = Multiply()([x, y])
    x_geno = np.ones([1, 1, y.shape[2]]) * geno


    # run convolutional layer on the output of the multiplication
    shared_conv = ConvBlock(79, 13, 4)
    x = shared_conv(x)
    x_geno = shared_conv(x_geno)
    print(x.shape)



    x = Flatten()(x)
    x_geno = Flatten()(x_geno)
    print("Flattened normal:", x.shape)
    print("geno:",x_geno.shape)



    shared_dense = layers.Dense(units=80, activation="relu",
                             # kernel_regularizer=regularizers.L2(L2)
                             )
    shared_drop = layers.Dropout(0.16582)

    x = shared_dense(x)
    x_geno = shared_dense(x_geno)
    x = shared_drop(x)
    x_geno = shared_drop(x_geno)


    pheno_out = Dense(1, activation="sigmoid", name="pheno_out", kernel_initializer=tf.keras.initializers.GlorotUniform())(x)

    gr = GradientReversalLayer()(x_geno)
    pheno_out_geno = Dense(1, activation="sigmoid", name="pheno_out_geno", kernel_initializer=tf.keras.initializers.GlorotUniform())(gr)

    pc_stop = Lambda(lambda x: K.stop_gradient(x))(x)
    pc_out = Dense(10, name="pc_out", activation="linear")(pc_stop)
    pc_to_pheno = Lambda(lambda x: K.stop_gradient(x))(pc_out)

    dense_layers = hp.Int("dense_layers", 1, 3)
    for k in range(dense_layers):
        with hp.conditional_scope("dense_layers", list(range(k + 1, dense_layers + 1))):
            units = hp.Int("dense_" + str(k), 5, 100)
            rate = hp.Float("dense_" + str(k) + "_dropout", min_value=0.0, max_value=0.5)
            pc_to_pheno = layers.Dense(units=units, activation="relu",
                                        # kernel_regularizer=regularizers.L2(L2)
                                        )(pc_to_pheno)
            pc_to_pheno = layers.Dropout(rate)(pc_to_pheno)


    pc_pheno_out = Dense(1, activation="sigmoid", name="pc_pheno_out",
                         kernel_initializer=tf.keras.initializers.GlorotUniform())(pc_to_pheno)


    losses = {
        "pheno_out": "binary_crossentropy",
        "pheno_out_geno": "binary_crossentropy",
        "pc_out" : "mse",
        "pc_pheno_out": "binary_crossentropy"
    }

    lossWeights = {
        "pheno_out": 1.0,
        "pheno_out_geno": 0.463,
        "pc_out" : 1,
        "pc_pheno_out": 1
    }

    model = Model(inputs=input, outputs=[pheno_out,
                                         pheno_out_geno,
                                         pc_out,
                                         pc_pheno_out
                  ])


    lr = 2.7617e-07

    model.compile(optimizer=keras.optimizers.Adam(learning_rate=lr),
                  loss=losses, loss_weights=lossWeights, metrics= {"pheno_out": "AUC", "pheno_out_geno": "AUC", "pc_out": tfa.metrics.RSquare(), "pc_pheno_out": "AUC"})



    return (model)


tuner = Hyperband(
    build_model,
    objective=keras_tuner.Objective("val_pc_pheno_out_auc_2", direction="max"),
    # max_trials=trials,
    # executions_per_trial=repeats,
    project_name="kt.adv_geno_stop_pc.11.11.22",
    hyperband_iterations= 2,
    max_epochs = epoch,
    # max_model_size = 100000,
    overwrite=True,
)
# early_stop = tf.keras.callbacks.EarlyStopping(monitor='val_loss', mode="min", patience=10)

# tuner.search_space_summary()

tuner.search(x=x_train, y ={"pheno_out": pheno_train, "pheno_out_geno": pheno_train, "pc_out": pca_train, "pc_pheno_out": pheno_train },
             validation_data=(x_valid, {"pheno_out": pheno_valid, "pheno_out_geno": pheno_valid, "pc_out": pca_valid, "pc_pheno_out": pheno_valid}),
             validation_batch_size=batch,
             validation_steps=(valid_size // batch),
             epochs=epoch,
             batch_size=batch,
             steps_per_epoch=(train_size // batch),
             )

tuner.results_summary()

best_model = tuner.get_best_models(num_models=1)[0]
best_model.summary()
