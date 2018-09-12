# Do-It-Yourself Single Camera 3D Pointer Input Device
**Started Spring 2016**

**University of Alberta, Department of Computing Science**

## Publication
B. Llanos, and Y.-H. Yang, “Do-It-Yourself Single Camera 3D Pointer Input Device,”
  in 2018 15th Conference on Computer and Robot Vision (CRV), 2018, pp. 214–221.
  DOI 10.1109/CRV.2018.00038

## Code Contributors
- Bernard Llanos, supervised by Dr. Y.-H. Yang

## Other Project Contributors
- CMPUT615 Project partners (Winter 2017 offering):
  - Kevin Gordon
  - Noah Weninger

## Description

Vision-based localization of a thin probe:
- Creation of a colour-based probe detection model.
- Detection of the probe in an image
- Determination of the probe's 3D position from a single camera's view
- Quantitative evaluation of probe localization

## Environment
- MATLAB R2016b and later were used for development, but earlier releases
  will likely work.
- Requires the Image Processing Toolbox, the Computer Vision Systems Toolbox,
  and the Statistics and Machine Learning toolbox (just for PCA).
- Live video processing uses the MATLAB Support Package for USB Webcams,
  and the Image Acquisition Toolbox in combination with the
  Image Acquisition Toolbox Support Package for Point Grey Hardware,
  depending on the type of camera to be used.

## Usage Instructions
- Files starting with capital letters are MATLAB scripts, which can be used
  as described in their documentation comments.
  - The scripts in 'data_analysis/' require data not provided in this repository,
    and most of them have no use beyond producing some of the figures in the
    conference paper cited above.
- Remaining files are MATLAB functions called by the scripts.

================================================================================

# Step-by-Step Walkthrough

1. Create a probe from a thin object by adding an arrangement of coloured bands.
   Refer to 'CreateProbeDetectionModel.m' for design details.

2. Measure the physical lengths and widths of the probe, as explained in
   'CreateProbeDetectionModel.m', and provide labels for its colours.

3. Calibrate a camera to determine its intrinsic and extrinsic parameters, and
   capture images of the probe with the camera.

4. Add this directory and all subdirectories to the MATLAB path.

5. Select one or more images, or capture a video, to calibrate the camera's
   RGB noise. Either use images or a video of a stationary scene as input
   for the script 'noise_estimation/EstimateRGBStandardDeviationsNonInteractive.m',
   or use arbitrary images as input for the script
   'noise_estimation/EstimateRGBStandardDeviationsInteractive.m'

6. Select an image, showing the entire probe, as the template for probe detection.
   Annotate the points where the coloured bands meet the probe contour in the image.
   Also mark the tip and the other end (which may or may not taper to a point).
   Refer to 'CreateProbeDetectionModel.m' for details.

7. Run 'CreateProbeDetectionModel.m' and save its output for use by the detector.

## Single-Image Workflow

8. Correct all images for lens distortion.

9. Choose an image in which to detect the probe.

10. Run 'DetectProbeScript.m' with the image as one of the inputs.
   If detection is unsuccessful, provide a bounding box by setting
   `request_bounding_region` to `true` and re-running the script.

11. Run 'LocalizeProbeScript.m' on the output of 'DetectProbeScript.m' to
    find the 3D pose of the probe.

## Live or Pre-Recorded Video Workflow

8. Use the script 'ProcessVideos.m'. Output CSV files (containing point clouds)
   can be visualized using the point cloud viewer at
   https://github.com/bllanos/point-cloud-viewer
