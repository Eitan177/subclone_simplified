import numpy as np
import pandas as pd
from Bio import pairwise2
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import streamlit as st


def framework_merge(framework_1, framework_2):
    """
    Merges two files between the same framework. Each index in framework number 1 is getting assigned correct merges from framework number 2
    :param framework_1: Arbitrary framework 1
    :param framework_2: Arbitrary framework 2
    :return: correct_merge: A list containing lists of matched merge sequences from framework 2 to each individual framework 1 sequence.
    """

    correct_merge = []
    framework_1_outliers = framework_1['Outliers'].tolist()
    framework_2_outliers = framework_2['Outliers'].tolist()

    for outlier1_position in range(len(framework_1_outliers)):
        position_merge = []
        for outlier2_position in range(len(framework_2_outliers)):
            shortest_length = min(len(framework_1_outliers[outlier1_position]), len(framework_2_outliers[outlier2_position]))
            alignments = pairwise2.align.globalms(framework_1_outliers[outlier1_position].replace("-", ""), framework_2_outliers[outlier2_position].replace("-", ""), 1, -1, -1, 0, score_only=True)
            if shortest_length - alignments <= 2:
                position_merge.append(outlier2_position)
        correct_merge.append(position_merge)
    print(correct_merge)
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


def heatmap_matrix_builder(outliers1, outliers2=None):
    """
    Builds the matrix for the heatmap matrix. Input either takes in a set of 3 frameworks containing outliers and does an
    interallele merging and creates a heatmap or input takes 2 sets of 3 frameworks of outliers and does a replicate
    comparison on them.
    :param outliers1: Outliers of 3 frameworks combined.
    :param outliers2: Optional. Outliers from 3
    :return:
    """

    if outliers2 is None:
        outliers2 = outliers1
    outlier_length_1 = len(outliers1)
    outlier_length_2 = len(outliers2)
    matrix = np.zeros((outlier_length_1, outlier_length_2))
    for outlier1_position in range(outlier_length_1):
        for outlier2_position in range(outlier_length_2):
            shortest_length = min(len(outliers1[outlier1_position].replace("-", "")), len(outliers2[outlier2_position].replace("-", "")))
            alignments = pairwise2.align.globalms(outliers1[outlier1_position].replace("-", ""), outliers2[outlier2_position].replace("-", ""), 1, -1, -1, 0, score_only=True)
            matrix[outlier1_position][outlier2_position] = shortest_length - alignments
    return matrix


def build_heatmap_graph(heatmap_matrix):
    fig, ax = plt.subplots()
    sns.heatmap(heatmap_matrix, linewidth=0.5)
    plt.title("Pairwise alignment score difference between two sequences")
    plt.xlabel("Sequence Indexes")
    plt.ylabel("Sequence Indexes")
    st.write(fig)
    sns.set()
    matrix_merge_threshold = np.where(heatmap_matrix >= 2, 1, 0)
    fig1, ax1 = plt.subplots()
    sns.heatmap(matrix_merge_threshold, linewidth=0.5)
    plt.title("Merges between two sequences")
    plt.xlabel("Sequence Indexes")
    plt.ylabel("Sequence Indexes")
    st.write(fig1)
    sns.set()


def Kevinsmerge(outliers_inter_allele_1, outliers_inter_allele_2):
    matrix = heatmap_matrix_builder(outliers_inter_allele_1, outliers_inter_allele_2)
    build_heatmap_graph(matrix)
