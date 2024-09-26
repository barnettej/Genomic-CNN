import os
import tensorflow as tf
import math
from tensorflow import keras
import numpy as np
from tensorflow import data
import keras_tuner
import pandas as pd
from tensorflow.keras import Input, optimizers, regularizers
from tensorflow.keras.layers import Dense, concatenate, Lambda, Multiply, Conv2D, Conv1D, Add, MaxPool1D, Dropout, Reshape, Flatten, \
    MaxPool2D, BatchNormalization, GlobalAveragePooling2D, GlobalAveragePooling1D
from tensorflow.keras import Model
from tensorflow.keras import layers
import keras.backend as K
import datetime
from tensorflow.keras.utils import plot_model
import itertools
from keras_tuner import BayesianOptimization, RandomSearch, Hyperband


batch = 
epoch = 
train_size = 
valid_size = 
test_size = 
trials = 
repeats = 

initializer = tf.keras.initializers.HeUniform()

train_path = "/home/user/user/ukbiobank/t2d/data/train.pca.csv"
valid_path = "/home/user/user/ukbiobank/t2d/data/valid.pca.csv"

train_ds = tf.data.experimental.make_csv_dataset(
    train_path,
    batch_size=batch,
    label_name="t2d",
    # num_epochs = 1,
    shuffle=True,
    num_parallel_reads=100,
    shuffle_buffer_size = 100,
    # shuffle_buffer_size=train_size,
    # select_columns = select_columns,
    ignore_errors=True)

valid_ds = tf.data.experimental.make_csv_dataset(
    valid_path,
    batch_size=batch,
    label_name="t2d",
    # num_epochs = 1,
    shuffle=True,
    num_parallel_reads=100,
    shuffle_buffer_size = 100,
    # shuffle_buffer_size=valid_size,
    # select_columns = select_columns,
    ignore_errors=True)



for feature_batch, label_batch in test_ds.take(1):
    # print("'hypertension': {}".format(label_batch))
    # print("features:")
    features = []
    for key, value in feature_batch.items():
        features.append(key)


def combine_batch_samples(samples, targets):
    inp1 = []
    for k in features:
        inp1.append(samples[k])

    # print(inp1)
    inp1 = tf.stack(inp1, axis=-1)

    return {'input_1': inp1}, targets


train_dataset = train_ds.map(combine_batch_samples)
train_dataset = train_dataset.take(train_size)
valid_dataset = valid_ds.map(combine_batch_samples)
valid_dataset = valid_dataset.take(valid_size)
test_dataset = test_ds.map(combine_batch_samples)
test_dataset = test_dataset.take(test_size)


def build_model(hp):
    dense_layers = hp.Int("dense_layers", 1, 3)

    # input genotype data
    input = Input(shape=(len(features)), batch_size=len(label_batch))
    x = input

    L2_bool = hp.Boolean("L2_bool")
    L2 = 0
    with hp.conditional_scope("L2_bool", ["True"]):
        if L2_bool == True:
            L2 = hp.Float("L2", min_value=1e-8, max_value=1e-6, sampling="log")

    for k in range(dense_layers):
        with hp.conditional_scope("dense_layers", list(range(k + 1, dense_layers + 1))):
            units = hp.Int("dense_" + str(k), 5, 100)
            rate = hp.Float("dense_" + str(k) + "_dropout", min_value=0.0, max_value=0.7)
            x = layers.Dense(units=units, activation="relu", kernel_regularizer=regularizers.L2(L2))(x)
            x = layers.Dropout(rate)(x)

    ht_out = Dense(1, activation="sigmoid", name="ht_out", kernel_initializer=tf.keras.initializers.GlorotUniform())(x)

    losses = {
        "ht_out": "binary_crossentropy",
    }

    lossWeights = {
        "ht_out": 1.0}

    model = Model(inputs=input, outputs=[ht_out], )

    lr = hp.Float("lr", min_value=1e-6, max_value=1e-3, sampling="log")

    model.compile(optimizer=keras.optimizers.Adam(learning_rate=lr),
                  loss=losses, loss_weights=lossWeights, metrics=["AUC", "accuracy"])


    return (model)


tuner = Hyperband(
    build_model,
    objective=keras_tuner.Objective("val_auc", direction="max"),
    # max_trials=trials,
    # executions_per_trial=repeats,
    project_name="t2d.nn.pca.11.10.22",
    hyperband_iterations= 5,
    # max_epochs = epoch,
    # max_model_size = 500000,
    overwrite=False,
)
# early_stop = tf.keras.callbacks.EarlyStopping(monitor='val_loss', mode="min", patience=5)

# tuner.search_space_summary()


tuner.search(train_dataset,
             validation_data=valid_dataset,
             validation_batch_size=batch,
             validation_steps=(valid_size // batch),
             epochs=epoch,
             batch_size=batch,
             steps_per_epoch=(train_size // batch),
             )

tuner.results_summary()

best_model = tuner.get_best_models(num_models=1)[0]
best_model.summary()
