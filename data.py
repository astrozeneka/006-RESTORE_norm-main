import numpy as np
from collections import defaultdict

def nested_dict():
    """
    A nested dictionary for hierarchical storage of thresholds.
    """
    return defaultdict(nested_dict)


def process_data(data, pos, neg, floor, pct_x=0.75, pct_y=0.75, sigma_x=3, sigma_y=3):
    x_outliers = data[pos] > data[pos].quantile(0.95) #(np.median(data[pos]) + sigma_x * data[pos].std())
    y_outliers = data[neg] > data[neg].quantile(0.95) #(np.median(data[neg]) + sigma_y * data[neg].std())
    joint_outliers = np.invert(x_outliers & y_outliers)

    x_inliers = data[pos] > np.quantile(data[pos], pct_x)
    y_inliers = data[neg] > np.quantile(data[neg], pct_y)
    joint_inliers = x_inliers | y_inliers
    
    x_floor = data[pos] > floor
    y_floor = data[neg] > floor
    joint_floor = x_floor & y_floor

    return data[joint_outliers & joint_inliers & joint_floor]
