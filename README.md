# Matlab_door_detection:

This Matlab script is designed to detect and measure door gaps in images using image processing techniques. It includes image preprocessing, feature extraction, and Hough transform for line detection.

Table of Contents:

    Introduction
    Requirements
    Usage
    Tasks Overview
    Example Test Images with Solutions

Introduction:

The script processes images to detect door contours and gaps by applying several image processing techniques. It includes steps such as contrast adjustment, noise reduction, gradient calculation, and Hough transform to identify vertical and horizontal lines representing door edges.
Requirements

    Matlab with Image Processing Toolbox

Usage:

    Clone the repository and navigate to the script directory.
    Place your test images in the same directory or adjust the image path in the script.
    Run the script in Matlab:

    MatLab_door_detection.m

Tasks Overview:

Task 1: Image Preprocessing

    Conversion to Gray Value Image: Converts the input image to grayscale.
    Adjust Contrast: Enhances the image contrast.
    Noise Reduction: Applies binomial low-pass filtering to reduce noise.

Task 2: Feature Extraction

    Gradient Calculation: Uses optimized Sobel operators to compute gradients in x and y directions.
    Edge Strength: Identifies strong edges by thresholding the gradient magnitudes.

Task 3: Segmentation of Door Gap Contour

    Vertical and Horizontal Edge Detection: Combines positive and negative gradients and applies morphological operations to detect continuous edges.

Task 4: Measuring Lines of the Door Gap

    Hough Transform: Computes the Hough transform for detected edges to identify and plot vertical and horizontal lines representing door edges.

Task 5: Measure and Classify Corner Points

    Corner Detection: Identifies intersection points of the detected lines to locate door corners.

Example Images with Solutions:

The script includes example images for testing:

    01 - R2441 - i.JPG
    01 - R2442 - i.JPG

Ensure your images are named and placed correctly in the directory.

To visualize the results, the script displays various figures showing intermediate and final processing steps.

Solutions:

![test_image_door_1](https://github.com/Paradox3333/Matlab_door_detection/assets/134848500/be23cef9-f31d-497b-9425-624502b01f35)



![test_image_door_2](https://github.com/Paradox3333/Matlab_door_detection/assets/134848500/d8106210-28c5-4869-a313-3db5a30c6919)
