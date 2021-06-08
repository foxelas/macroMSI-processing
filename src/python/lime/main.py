# -*- coding: utf-8 -*-

import tensorflow as tf
slim = tf.contrib.slim
import sys
sys.path.append('D:/elena/models/tf-models/slim')
#%load_ext autoreload
#%autoreload 2
import matplotlib.pyplot as plt
#%matplotlib inline
import numpy as np

from nets import inception
from preprocessing import inception_preprocessing
session = tf.Session()
image_size = inception.inception_v3.default_image_size

def transform_img_fn(path_list):
    out = []
    for f in path_list:
        image_raw = tf.image.decode_jpeg(open(f).read(), channels=3, encoding='cp932')
        image = inception_preprocessing.preprocess_image(image_raw, image_size, image_size, is_training=False)
        out.append(image)
    return session.run([out])[0]

from datasets import imagenet
names = imagenet.create_readable_names_for_imagenet_labels()
processed_images = tf.placeholder(tf.float32, shape=(None, 299, 299, 3))
import os
with slim.arg_scope(inception.inception_v3_arg_scope()):
    logits, _ = inception.inception_v3(processed_images, num_classes=1001, is_training=False)
probabilities = tf.nn.softmax(logits)

checkpoints_dir = 'D:/elena/models/tf-models/slim/pretrained'
init_fn = slim.assign_from_checkpoint_fn(
    os.path.join(checkpoints_dir, 'inception_v3.ckpt'),
    slim.get_model_variables('InceptionV3'))
init_fn(session)
def predict_fn(images):
    return session.run(probabilities, feed_dict={processed_images: images})

images = transform_img_fn(['D:/elena/Google Drive/titech/research/experiments/input/py_tests/dogs.jpg'])
# I'm dividing by 2 and adding 0.5 because of how this Inception represents images
plt.imshow(images[0] / 2 + 0.5)
preds = predict_fn(images)
for x in preds.argsort()[0][-5:]:
    print (x, names[x], preds[0,x])