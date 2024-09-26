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
from sklearn.preprocessing import LabelEncoder

batch = 
epoch = 
train_size = 
valid_size = 
test_size = 
trials = 
repeats = 

initializer = tf.keras.initializers.HeUniform()

train_path = "/home/user/user/ukbiobank/t2d/data/t2d.annocut.p0.01.train.pca.scale.csv"
valid_path = "/home/user/user/ukbiobank/t2d/data/t2d.annocut.p0.01.valid.pca.scale.csv"

train_ds = pd.read_csv(train_path)
print(train_ds.head())
x_train, pheno_train = train_ds.values[:55168,:10], train_ds.values[:55168,-1]

x_train = tf.cast(tf.convert_to_tensor(x_train), tf.float32)
pheno_train = LabelEncoder().fit_transform(pheno_train)

valid_ds = pd.read_csv(valid_path)
x_valid, pheno_valid = valid_ds.values[:11712,:10], valid_ds.values[:11712,-1]
x_valid = tf.cast(tf.convert_to_tensor(x_valid), tf.float32)
pheno_valid = LabelEncoder().fit_transform(pheno_valid)



print("x shape:",x_train.shape)
print("pheno shape:", pheno_train.shape)


def build_model(hp):
    dense_layers = hp.Int("dense_layers", 1, 3)

    # input genotype data
    input = Input(shape=(x_train.shape[1]), batch_size=batch)
    x = input


    for k in range(dense_layers):
        with hp.conditional_scope("dense_layers", list(range(k + 1, dense_layers + 1))):
            units = hp.Int("dense_" + str(k), 5, 100)
            rate = hp.Float("dense_" + str(k) + "_dropout", min_value=0.0, max_value=0.7)
            x = layers.Dense(units=units, activation="relu",
                             # kernel_regularizer=regularizers.L2(L2)
                             )(x)
            x = layers.Dropout(rate)(x)

    pheno_out = Dense(1, activation="sigmoid", name="pheno_out", kernel_initializer=tf.keras.initializers.GlorotUniform())(x)

    losses = {
        "pheno_out": "binary_crossentropy",
    }

    lossWeights = {
        "pheno_out": 1.0,
    }

    model = Model(inputs=input, outputs=[pheno_out], )

    lr = hp.Float("lr", min_value=1e-6, max_value=1e-3, sampling="log")


    model.compile(optimizer=keras.optimizers.Adam(learning_rate=lr),
                  loss=losses, loss_weights=lossWeights, metrics= {"pheno_out": "AUC"})


    return (model)


tuner = Hyperband(
    build_model,
    objective=keras_tuner.Objective("val_auc", direction="max"),
    # max_trials=trials,
    # executions_per_trial=repeats,
    project_name="kt.nn.pca_pred_pheno.11.10.22",
    hyperband_iterations= 3,
    # max_epochs = epoch,
    # max_model_size = 500000,
    overwrite=True,
)
# early_stop = tf.keras.callbacks.EarlyStopping(monitor='val_loss', mode="min", patience=5)

# tuner.search_space_summary()


tuner.search(x=x_train, y ={"pheno_out": pheno_train},
             validation_data=(x_valid, {"pheno_out": pheno_valid}),
             validation_batch_size=batch,
             validation_steps=(valid_size // batch),
             epochs=epoch,
             batch_size=batch,
             steps_per_epoch=(train_size // batch),
             )

tuner.results_summary()

best_model = tuner.get_best_models(num_models=1)[0]
best_model.summary()
