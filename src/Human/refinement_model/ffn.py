from tensorflow.keras import backend as K
import tensorflow as tf
from tensorflow.keras import layers
import numpy as np
from tensorflow.keras.utils import plot_model
from tensorflow.keras.models import load_model
from tensorflow.keras.callbacks import ModelCheckpoint
import matplotlib.pyplot as plt

class FFN:
    def __init__(self, input_shape, train_feature, train_label,
                 epochs=2000, s=None,
                 eval_=True, learning_curve=False):
        self.input_shape = input_shape
        self.train_feature = train_feature
        self.train_label = train_label
        self.epochs = epochs
        self.s = s
        self.eval_ = eval_
        self.learning_curve = learning_curve
        # self.fold = fold

    def run(self):
        model = self.init_model()
        model.fit(self.train_feature, self.train_label,
                  batch_size=np.size(a=self.train_feature, axis=0),
                  epochs=self.epochs, verbose=0)
        return model

    def init_model(self):
        # 1) Initialize Sequential
        model = tf.keras.Sequential()

        # 2) Add layer
        model.add(layers.Dense(units=16, input_dim=self.input_shape))
        model.add(layers.BatchNormalization())
        model.add(layers.Dropout(0.2))
        model.add(layers.ReLU())

        model.add(layers.Dense(32, kernel_regularizer=tf.keras.regularizers.l2(0.01)))
        model.add(layers.BatchNormalization())
        model.add(layers.Dropout(0.2))
        model.add(layers.ReLU())

        model.add(layers.Dense(128, kernel_regularizer=tf.keras.regularizers.l2(0.01)))
        model.add(layers.BatchNormalization())
        model.add(layers.Dropout(0.2))
        model.add(layers.ReLU())

        model.add(layers.Dense(64, kernel_regularizer=tf.keras.regularizers.l2(0.01)))
        model.add(layers.BatchNormalization())
        model.add(layers.Dropout(0.2))
        model.add(layers.ReLU())

        model.add(layers.Dense(32, kernel_regularizer=tf.keras.regularizers.l2(0.01)))
        model.add(layers.BatchNormalization())
        model.add(layers.Dropout(0.2))
        model.add(layers.ReLU())

        model.add(layers.Dense(units=1, activation='sigmoid', kernel_regularizer=tf.keras.regularizers.l2(0.01)))

        # 3) Initialize cost function, optimization algorithm, and metrics
        model.compile(optimizer='adam', loss='binary_crossentropy')

        model.epochs = self.epochs

        return model

    def accuracy(self, target, pred):
        target = K.round(K.clip(x=target, min_value=0, max_value=1))
        pred = K.round(K.clip(pred, 0, 1))
        accuracy = K.mean(K.equal(target, pred))
        return accuracy

    def recall(self, target, pred):
        """ Calculate recall
            Args:
                target: target values of data
                pred: prediction values of model
            1) Clip target and pred appropriately
            2) Calculate true positive
            3) Calculate false negative + true positive
            4) Calculate recall
        """
        # 1) Clip target and pred appropriately
        target = K.round(K.clip(x=target, min_value=0, max_value=1))
        pred = K.round(K.clip(pred, 0, 1))
        """ tf.keras.backend.clip
                Element-wise value clipping.
                Args:
                    x: Tensor or variable.
                    min_value: Python float or integer.
                    max_value: Python float or integer.
        """
        """ tf.keras.backend.round
                Element-wise rounding to the closest integer.
                In case of tie, the rounding mode used is "half to even".
                Args:
                    x: Tensor or variable.
        """
        # 2) Calculate true positive
        # Recall that every value is 1 or 0 so, sum means amount of data which is 1.
        # And target * pred means AND(target, pred) so, they activate only when they are true positive
        n_true_positive = K.sum(target * pred)

        # 3) Calculate false negative + true positive
        n_true_positive_false_negative = K.sum(target)

        # 4) Calculate recall
        # recall =  (true Positive) / (true Positive + false Negative)
        # We add very small value by using K.epsilon() to prevent division by zero error
        recall = n_true_positive / (n_true_positive_false_negative + K.epsilon())
        return recall

    def precision(self, target, pred):
        """ Calculate precision
            Args:
                target: target values of data
                pred: prediction values of model
            1) Clip target and pred appropriately
            2) Calculate true positive
            3) Calculate amount of retrieved data
            4) Calculate precision
        """
        # 1) Clip target and pred appropriately
        target = K.round(K.clip(target, 0, 1))
        pred = K.round(K.clip(pred, 0, 1))
        """ tf.keras.backend.clip
                    Element-wise value clipping.
                    Args:
                        x: Tensor or variable.
                        min_value: Python float or integer.
                        max_value: Python float or integer.
            """
        """ tf.keras.backend.round
                Element-wise rounding to the closest integer.
                In case of tie, the rounding mode used is "half to even".
                Args:
                    x: Tensor or variable.
        """
        # 2) Calculate true positive
        # Recall that every value is 1 or 0 so, sum means amount of data which is 1.
        # And target * pred means AND(target, pred) so, they activate only when they are true positive
        n_true_positive = K.sum(target * pred)

        # 3) Calculate amount of retrieved data
        # amount of retrieved data = false positive + true positive
        n_retrieved_data = K.sum(pred)

        # 4) Calculate precision
        # precision = (true Positive) / (false positive + true positive)
        # We add very small value by using K.epsilon() to prevent division by zero error
        precision = n_true_positive / (n_retrieved_data + K.epsilon())
        return precision

    def f1(self, target, pred):
        precision = self.precision(target, pred)
        recall = self.precision(target, pred)

        f1 = 2 * (recall * precision) / (recall + precision)
        return f1
