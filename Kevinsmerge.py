import numpy as np
import pandas as pd
from Bio import pairwise2


def framework_merge(framework_1, framework_2):
    """
    Merges two files between the same framework. Each index in framework number 1 is getting assigned correct merges from framework number 2
    :param framework_1: Arbitrary framework 1
    :param framework_2: Arbitrary framework 2
    :return: correct_merge: A list containing lists of matched merge sequences from framework 2.
    """

    correct_merge = []
    framework_1_outliers = framework_1['seq'].tolist()
    framework_2_outliers = framework_2['seq'].tolist()

    for outlier1_position in range(len(framework_1_outliers)):
        position_merge = []
        for outlier2_position in range(len(framework_2_outliers)):
            shortest_length = min(len(framework_1_outliers[outlier1_position]), len(framework_2_outliers[outlier2_position]))
            alignments = pairwise2.align.globalms(framework_1_outliers[outlier1_position].replace("-", ""), framework_2_outliers[outlier2_position].replace("-", ""), 1, -1, -1, 0, score_only=True)
            if shortest_length - alignments <= 2:
                position_merge.append(outlier2_position)
        correct_merge.append(position_merge)
    return correct_merge


def merge_outliers(outliers):
    """
    Appends 2 columns comparing against the 2 other framework for matches. Each column contains a list of indexes from the outlier file
    :param outliers: The dataframe of the 3 frameworks.
    :return: The dataframe of each framework with appended column comparing against index of other framework for matches.
             Column contains list of indexes from other framework that match.
    """

    merged_outliers = []
    for framework_number_1 in range(3):
        merges = []
        merge_column_name = []
        for framework_number_2 in range(3):
            if framework_number_1 != framework_number_2:
                merge_column_name.append('framework_' + str(framework_number_1) + '_merge_with_framework_' + str(framework_number_2))
                merges.append(framework_merge(outliers[framework_number_1], outliers[framework_number_2]))
        outliers[framework_number_1][merge_column_name[0]] = merges[0]
        outliers[framework_number_1][merge_column_name[1]] = merges[1]
        merged_outliers.append(outliers[framework_number_1])
    return merged_outliers


def merge(framework1,framework2,framework3):

    outliers = [framework1, framework2, framework3]
    framework1_merge, framework2_merge, framework3_merge = merge_outliers(outliers)
    framework1_merge.to_csv("fr1_merge.csv", index=False)
    framework2_merge.to_csv("f2_merge.csv", index=False)
    framework3_merge.to_csv("fr3_merge.csv", index=False)
    return [framework1_merge, framework2_merge, framework3_merge]
