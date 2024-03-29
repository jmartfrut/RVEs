
function visualizeMesh()
    # Launch the GUI to see the results:
    if !("-nopopup" in ARGS)
        gmsh.fltk.run()
    end

end

function stopGmsh()
    Gmsh.finalize()
end

function startGmsh()
    gmsh.initialize()
end

function ShowInfo(a::Int64)
    gmsh.option.setNumber("General.Terminal",a)
end

function saveMesh(output_file::String)
    gmsh.write(output_file)
end