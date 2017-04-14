# CMPUT615 Project and Ongoing Research Endeavour
**Started Spring 2016**
**University of Alberta, Department of Computing Science**

## Contributors
- Bernard Llanos, supervised by Dr. Y.H. Yang
- CMPUT615 Project partners (Winter 2017 offering):
  - Kevin Gordon
  - Noah Weninger

## Description

Vision-based localization of a thin probe:
- Creation of a colour-based probe detection model.
- Detection of the probe in an image
- Determination of the probe's 3D position from a single camera's view
- Quantitative evaluation of probe localization

## Usage Instructions
- Files starting with capital letters are MATLAB scripts, which can be used
  as described in their documentation comments.
- Remaining files are MATLAB functions called by the scripts.

================================================================================

# Step-by-Step Walkthrough

1. Create a probe from a thin object by adding an arrangement of coloured bands.
   Refer to 'CreateProbeDetectionModel.m' for design details.

2. Measure the physical lengths and widths of the probe, as explained in
   'CreateProbeDetectionModel.m', and provide labels for its colours.

3. Calibrate a camera to determine its intrinsic and extrinsic parameters, and
   capture images of the probe with the camera.

4. Select one or more images to calibrate the camera's RGB noise, and use them
   as input for the script 'EstimateRGBStandardDeviations.m'.

5. Correct all images for lens distortion.

6. Select an image, showing the entire probe, as the template for probe detection.
   Annotate the points where the coloured bands meet the probe contour in the image.
   Also mark the tip and the other end (which may or may not taper to a point).
   Refer to 'CreateProbeDetectionModel.m' for details.

7. Run 'CreateProbeDetectionModel.m' and save its output for use by the detector.

8. Choose an image in which to detect the probe.

9. Run 'DetectProbe.m' with the image as one of the inputs. If detection is
   unsuccessful, provide a bounding box by setting `request_bounding_region` to `true`.

10. Run 'LocalizeProbe.m' on the output of 'DetectProbe.m' to find the 3D pose
    of the probe.
