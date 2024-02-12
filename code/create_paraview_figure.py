# trace generated using paraview version 5.11.0
# import paraview
# paraview.compatibility.major = 5
# paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *


def main(filename, min_value, max_value, name, label, blocknr, figure_folder):
    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'Xdmf3ReaderT'
    resultsxdmf = Xdmf3ReaderT(
        registrationName="results.xdmf",
        FileName=[
            filename,
        ],
    )

    # get animation scene
    animationScene1 = GetAnimationScene()

    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()

    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # get active view
    renderView1 = GetActiveViewOrCreate("RenderView")

    # # create a new 'Extract Block'
    extractBlock = ExtractBlock(registrationName="ExtractBlock", Input=resultsxdmf)

    # # Properties modified on extractBlock
    extractBlock.Selectors = [f"/Root/Block{blocknr}"]

    # show data in view
    extractBlockDisplay = Show(
        extractBlock, renderView1, "UnstructuredGridRepresentation"
    )

    # get color transfer function/color map for 'sigma_xx'
    sigma_xxLUT = GetColorTransferFunction(name)

    # get opacity transfer function/opacity map for 'sigma_xx'
    sigma_xxPWF = GetOpacityTransferFunction(name)

    # trace defaults for the display properties.
    extractBlockDisplay.Representation = "Surface"
    extractBlockDisplay.ColorArrayName = ["POINTS", name]
    extractBlockDisplay.LookupTable = sigma_xxLUT
    extractBlockDisplay.SelectTCoordArray = "None"
    extractBlockDisplay.SelectNormalArray = "None"
    extractBlockDisplay.SelectTangentArray = "None"
    extractBlockDisplay.OSPRayScaleArray = name
    extractBlockDisplay.OSPRayScaleFunction = "PiecewiseFunction"
    extractBlockDisplay.SelectOrientationVectors = "None"
    extractBlockDisplay.ScaleFactor = 400.0
    extractBlockDisplay.SelectScaleArray = name
    extractBlockDisplay.GlyphType = "Arrow"
    extractBlockDisplay.GlyphTableIndexArray = name
    extractBlockDisplay.GaussianRadius = 20.0
    extractBlockDisplay.SetScaleArray = ["POINTS", name]
    extractBlockDisplay.ScaleTransferFunction = "PiecewiseFunction"
    extractBlockDisplay.OpacityArray = ["POINTS", name]
    extractBlockDisplay.OpacityTransferFunction = "PiecewiseFunction"
    extractBlockDisplay.DataAxesGrid = "GridAxesRepresentation"
    extractBlockDisplay.PolarAxes = "PolarAxesRepresentation"
    extractBlockDisplay.ScalarOpacityFunction = sigma_xxPWF
    extractBlockDisplay.ScalarOpacityUnitDistance = 120.50872400061667
    extractBlockDisplay.OpacityArrayName = ["POINTS", name]
    extractBlockDisplay.SelectInputVectors = [None, ""]
    extractBlockDisplay.WriteLog = ""

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    extractBlockDisplay.ScaleTransferFunction.Points = [
        2.5313084643816214e-16,
        0.0,
        0.5,
        0.0,
        2.531850565467864e-16,
        1.0,
        0.5,
        0.0,
    ]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    extractBlockDisplay.OpacityTransferFunction.Points = [
        2.5313084643816214e-16,
        0.0,
        0.5,
        0.0,
        2.531850565467864e-16,
        1.0,
        0.5,
        0.0,
    ]

    # hide data in view
    Hide(resultsxdmf, renderView1)

    # show color bar/color legend
    extractBlockDisplay.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # get 2D transfer function for 'sigma_xx'
    sigma_xxTF2D = GetTransferFunction2D(name)

    animationScene1.GoToNext()

    # Rescale transfer function
    sigma_xxLUT.RescaleTransferFunction(-2.8939850330352783, 7.81556510925293)

    # Rescale transfer function
    sigma_xxPWF.RescaleTransferFunction(-2.8939850330352783, 7.81556510925293)

    # create a new 'Clip'
    clip1 = Clip(registrationName="Clip1", Input=extractBlock)
    clip1.ClipType = "Plane"
    clip1.HyperTreeGridClipper = "Plane"
    clip1.Scalars = ["POINTS", name]
    clip1.Value = 2.4607900381088257

    # show data in view
    clip1Display = Show(clip1, renderView1, "UnstructuredGridRepresentation")

    # trace defaults for the display properties.
    clip1Display.Representation = "Surface"
    clip1Display.ColorArrayName = ["POINTS", name]
    clip1Display.LookupTable = sigma_xxLUT
    clip1Display.SelectTCoordArray = "None"
    clip1Display.SelectNormalArray = "None"
    clip1Display.SelectTangentArray = "None"
    clip1Display.OSPRayScaleArray = name
    clip1Display.OSPRayScaleFunction = "PiecewiseFunction"
    clip1Display.SelectOrientationVectors = "None"
    clip1Display.ScaleFactor = 200.0
    clip1Display.SelectScaleArray = name
    clip1Display.GlyphType = "Arrow"
    clip1Display.GlyphTableIndexArray = name
    clip1Display.GaussianRadius = 10.0
    clip1Display.SetScaleArray = ["POINTS", name]
    clip1Display.ScaleTransferFunction = "PiecewiseFunction"
    clip1Display.OpacityArray = ["POINTS", name]
    clip1Display.OpacityTransferFunction = "PiecewiseFunction"
    clip1Display.DataAxesGrid = "GridAxesRepresentation"
    clip1Display.PolarAxes = "PolarAxesRepresentation"
    clip1Display.ScalarOpacityFunction = sigma_xxPWF
    clip1Display.ScalarOpacityUnitDistance = 80.35903137837035
    clip1Display.OpacityArrayName = ["POINTS", name]
    clip1Display.SelectInputVectors = [None, ""]
    clip1Display.WriteLog = ""

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    clip1Display.ScaleTransferFunction.Points = [
        -2.8062829971313477,
        0.0,
        0.5,
        0.0,
        7.81556510925293,
        1.0,
        0.5,
        0.0,
    ]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    clip1Display.OpacityTransferFunction.Points = [
        -2.8062829971313477,
        0.0,
        0.5,
        0.0,
        7.81556510925293,
        1.0,
        0.5,
        0.0,
    ]

    # hide data in view
    Hide(extractBlock, renderView1)

    # show color bar/color legend
    clip1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source
    SetActiveSource(extractBlock)

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=clip1.ClipType)

    # create a new 'Clip'
    clip2 = Clip(registrationName="Clip2", Input=extractBlock)
    clip2.ClipType = "Plane"
    clip2.HyperTreeGridClipper = "Plane"
    clip2.Scalars = ["POINTS", name]
    clip2.Value = 2.4607900381088257

    # Properties modified on clip2.ClipType
    clip2.ClipType.Normal = [0.0, 1.0, 0.0]

    # show data in view
    clip2Display = Show(clip2, renderView1, "UnstructuredGridRepresentation")

    # trace defaults for the display properties.
    clip2Display.Representation = "Surface"
    clip2Display.ColorArrayName = ["POINTS", name]
    clip2Display.LookupTable = sigma_xxLUT
    clip2Display.SelectTCoordArray = "None"
    clip2Display.SelectNormalArray = "None"
    clip2Display.SelectTangentArray = "None"
    clip2Display.OSPRayScaleArray = name
    clip2Display.OSPRayScaleFunction = "PiecewiseFunction"
    clip2Display.SelectOrientationVectors = "None"
    clip2Display.ScaleFactor = 400.0
    clip2Display.SelectScaleArray = name
    clip2Display.GlyphType = "Arrow"
    clip2Display.GlyphTableIndexArray = name
    clip2Display.GaussianRadius = 20.0
    clip2Display.SetScaleArray = ["POINTS", name]
    clip2Display.ScaleTransferFunction = "PiecewiseFunction"
    clip2Display.OpacityArray = ["POINTS", name]
    clip2Display.OpacityTransferFunction = "PiecewiseFunction"
    clip2Display.DataAxesGrid = "GridAxesRepresentation"
    clip2Display.PolarAxes = "PolarAxesRepresentation"
    clip2Display.ScalarOpacityFunction = sigma_xxPWF
    clip2Display.ScalarOpacityUnitDistance = 145.3794450275547
    clip2Display.OpacityArrayName = ["POINTS", name]
    clip2Display.SelectInputVectors = [None, ""]
    clip2Display.WriteLog = ""

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    clip2Display.ScaleTransferFunction.Points = [
        -2.363538980484009,
        0.0,
        0.5,
        0.0,
        6.50651741027832,
        1.0,
        0.5,
        0.0,
    ]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    clip2Display.OpacityTransferFunction.Points = [
        -2.363538980484009,
        0.0,
        0.5,
        0.0,
        6.50651741027832,
        1.0,
        0.5,
        0.0,
    ]

    # hide data in view
    Hide(extractBlock, renderView1)

    # show color bar/color legend
    clip2Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # get color legend/bar for sigma_xxLUT in view renderView1
    sigma_xxLUTColorBar = GetScalarBar(sigma_xxLUT, renderView1)

    # Properties modified on sigma_xxLUTColorBar
    sigma_xxLUTColorBar.Title = label
    sigma_xxLUTColorBar.RangeLabelFormat = "%-#6.1f"

    # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
    sigma_xxLUT.ApplyPreset("Viridis (matplotlib)", True)

    # change scalar bar placement
    sigma_xxLUTColorBar.Orientation = "Horizontal"
    sigma_xxLUTColorBar.WindowLocation = "Any Location"
    sigma_xxLUTColorBar.Position = [0.3250956464379949, 0.7651456310679611]
    sigma_xxLUTColorBar.ScalarBarLength = 0.33

    # Rescale transfer function
    sigma_xxLUT.RescaleTransferFunction(min_value, max_value)

    # Rescale transfer function
    sigma_xxPWF.RescaleTransferFunction(min_value, max_value)

    # Rescale 2D transfer function
    sigma_xxTF2D.RescaleTransferFunction(min_value, max_value, 0.0, 1.0)

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=clip2.ClipType)

    # get layout
    layout1 = GetLayout()

    # --------------------------------
    # saving layout sizes for layouts

    # layout/tab size in pixels
    layout1.SetSize(1516, 1030)

    # -----------------------------------
    # saving camera placements for views

    # current camera placement for renderView1
    renderView1.CameraPosition = [
        4527.975828521292,
        1143.9996458085711,
        1335.8353866420669,
    ]
    renderView1.CameraFocalPoint = [
        -501.3331417934615,
        -126.6625437863954,
        -147.90241305746002,
    ]
    renderView1.CameraViewUp = [
        -0.21900790347787114,
        0.9715903578242485,
        -0.08970905638326834,
    ]
    renderView1.CameraParallelScale = 2044.5048300260873

    # export view
    ExportView(f"{figure_folder}/{name}.svg", view=renderView1)

    # ================================================================
    # addendum: following script captures some of the application
    # state to faithfully reproduce the visualization during playback
    # ================================================================

    # --------------------------------------------
    # uncomment the following to render all views
    # RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--filename", type=str)
    parser.add_argument("--min_value", type=float)
    parser.add_argument("--max_value", type=float)
    parser.add_argument("--name", type=str)
    parser.add_argument("--label", type=str)
    parser.add_argument("--blocknr", type=int)
    parser.add_argument("--figure_folder", type=str)
    args = parser.parse_args()

    main(
        filename=args.filename,
        min_value=float(args.min_value),
        max_value=float(args.max_value),
        name=args.name,
        label=args.label,
        blocknr=int(args.blocknr),
        figure_folder=args.figure_folder,
    )
