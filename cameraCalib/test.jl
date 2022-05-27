using Images, ImageView, Gtk, Gtk.ShortNames
using Spinnaker

```
SpinEnumNode SpinFloatNode SpinBooleanNode SpinIntegerNode
CameraNodeMap CameraTLStreamNodeMap CameraTLDeviceNodeMap
```

"""
2カメラの設定
"""
function configSetup(camList, exposure = 6001.0, gain = 13.10)
    cam1 = camList[0]
    cam2 = camList[1]
    
    gammenset1 = gammaenable!(cam1,false)
    @assert gammenset1 == false
    gammenset2 = gammaenable!(cam2,false)
    @assert gammenset2 == false

    adcbits!(cam1, "Bit12")
    @assert adcbits(cam1) == "Bit12" 
    adcbits!(cam2, "Bit12")
    @assert adcbits(cam2) == "Bit12"

    Spinnaker.set!(Spinnaker.SpinBooleanNode(cam1,"AcquisitionFrameRateEnable"),false)
    Spinnaker.set!(Spinnaker.SpinBooleanNode(cam2,"AcquisitionFrameRateEnable"),false)


    pixelformat!(cam1, "Mono16")
    @assert pixelformat(cam1) == "Mono16"
    pixelformat!(cam2, "Mono16")
    @assert pixelformat(cam2) == "Mono16"

    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"ExposureAuto"),"Off")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam2,"ExposureAuto"),"Off")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"ExposureMode"),"Timed")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"ExposureMode"),"Timed")
    exposure!(cam1, exposure)
    exposure!(cam2, exposure)

    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"GainAuto"),"Off")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam2,"GainAuto"),"Off")
    gain!(cam1, gain)
    gain!(cam2, gain)

    triggermode!(cam1, "On")
    triggersource!(cam1, "Software")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"TriggerSelector"),"FrameStart")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"LineSelector"),"Line1")
    Spinnaker.line_mode(cam1, "Output")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"LineSelector"),"Line2")
    Spinnaker.v3_3_enable(cam1,true)
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam2,"TriggerOverlap"),"Off")

    triggermode!(cam2, "On")
    triggersource!(cam2, "Line3")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam2,"TriggerSelector"),"FrameStart")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam2,"TriggerOverlap"),"ReadOut")

    Spinnaker.set!(Spinnaker.SpinEnumNode(cam1,"AcquisitionMode"),"Continuous")
    Spinnaker.set!(Spinnaker.SpinEnumNode(cam2,"AcquisitionMode"),"Continuous")
end

function imgAcquisition(camList,state, stopbutton, canvases)
    cam1 = camList[0]
    cam2 = camList[1]

    start!(cam1)
    start!(cam2)
    signal_connect(x -> state = "stop", stopbutton, "clicked")
    while state == "running"
        # start!(cam1)
        # start!(cam2)
        trigger!(cam1)
        img1 = getimage(cam1)
        img2 = getimage(cam2)
        arr1 = CameraImage(img1,Float64, normalize = true)
        arr2 = CameraImage(img2,Float64, normalize = true)
        println("img1 maximum: ",maximum(arr1))
        println("img2 maximum: ",maximum(arr2))
        println("")
        arr1[1,1] = 0.99
        arr2[1,1] = 0.99
        imshow(canvases[1,1], arr1)
        imshow(canvases[1,2], arr2)
        Gtk.showall(gui["window"])
        sleep(0.01)
        # stop!(cam1)
        # stop!(cam2)
    end
    stop!(cam1)
    stop!(cam2)
    Spinnaker._release!(cam1)
    Spinnaker._release!(cam2)
end

# get camera list
if !(@isdefined camList)
    camList = CameraList()
    println(camList)
    configSetup(camList)
end

win = Window("Camera Controller")
v = Box(:v)
l1 = Label("")
b1 = Button("Start")
l2 = Label("")
b2 = Button("Stop")
l3 = Label("")
b3 = Button("Finish")
l4 = Label("")
push!(win, v)
push!(v,l1)
push!(v, b1)
push!(v,l2)
push!(v, b2)
push!(v,l3)
push!(v, b3)
push!(v,l4)

showall(win)

gui = imshow_gui((1500,300),(1,2))
canvases = gui["canvas"]

imgAcquisTask(camList, state, stopbutton, canvases) = @task begin; imgAcquisition(camList,state,stopbutton,canvases); end

function start(camList, stopbutton, canvases)
    state = "running"
    schedule(imgAcquisTask(camList, state, stopbutton, canvases))
end

signal_connect(x -> start(camList, b2, canvases), b1, "clicked")
signal_connect(x -> exit(), b3, "clicked")

if !isinteractive()
    c = Condition()
    signal_connect(win, :destroy) do widget
        notify(c)
    end
    wait(c)
end