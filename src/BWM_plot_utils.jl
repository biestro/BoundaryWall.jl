using Makie


"""
  getPlottableMesh(ms::SimpleMesh)

Converts a `Meshes.SimpleMesh` object into plottable data for `Makie.mesh`.
Plot using Makie.mesh(tcoords, tmatrix, colors = eachindex(tcoords)).

Useful since Meshes.viz! exports as png, and for colorbars
"""
function getPlottableMesh(ms::SimpleMesh)
  _mesh = ms
  topo = _mesh.topology
  dim = embeddim(_mesh)
  nvert = nvertices(_mesh)
  nelem = nelements(_mesh)
  verts = vertices(_mesh)
  topo = topology(_mesh)
  elems = elements(topo)
  coords = coords.(verts)
  tris4elem = map(elems) do elem
    I = indices(elem)
    [[I[1], I[i], I[i + 1]] for i in 2:(length(I) - 1)]
  end

  # flatten vector of triangles
  tris = [tri for tris in tris4elem for tri in tris]

  tcoords = coords
  tconnec = tris
  tmatrix = reduce(hcat, tconnec) |> transpose

  return tcoords, tmatrix
  ## in order to plot Makie.mesh(tcoords, tmatrix, colors = eachindex(tcoords))
end

"""
  meshesToPolygons(M::SimpleMesh)

**Use `getPlottableMesh(ms::SimpleMesh)` instead**. Converts a `Meshes.SimpleMesh` 
into a `Makie.Polygon`, used for plotting in an `.svg` format.
"""
function meshesToPolygons(M::SimpleMesh)
  _verts = (getindex.(M.vertices,1), getindex.(M.vertices, 2))
  _faces = map((elem)->CartesianIndex.([elem.indices...]'), M.topology.connec)
  # get x, y polygons
  _x_points = vec.(getindex.(Ref(_verts[1]), _faces))
  _y_points = vec.(getindex.(Ref(_verts[2]), _faces))

  _makie_points = [Makie.Point2f.(_x, _y) for (_x, _y) in zip(_x_points, _y_points)]

  return Makie.Polygon.(_makie_points)
end




theme = Theme(
  Axis = (
    xgridvisible=false, 
    ygridvisible=false, 
    xtickalign=1, 
    ytickalign=1,
    xminorticksvisible = true,
    yminorticksvisible = true,
    xminorticks=IntervalsBetween(9),
    yminorticks=IntervalsBetween(9),
    xminortickalign=1,
    yminortickalign=1,
    xticksmirrored=true,
    xminorticksmirrored=true,
    yticksmirrored=true,
    yminorticksmirrored=true,


    ), 
  Lines=(; 
  cycle=:linestyle),
  # palette = (color = to_colormap(:seaborn_colorblind), marker = []),
  Scatter=(;
    cycle=:marker,
    # strokewidth=2,
    # color=RGBAf(0.,0.,0.,0.),
    strokecolor=:black),
  ScatterLines=(;
    cycle=Cycle([:linestyle, :marker], covary = true),
    ),

    )


# colors = to_colormap(:okabe_ito)
colors = Makie.wong_colors()
cmap = Makie.to_colormap(Reverse(:linear_protanopic_deuteranopic_kbw_5_98_c40_n256))
cmap[1] = RGBAf(1,1,1,1)
# cmap = Makie.to_colormap(Reverse(:linear_kbgyw_5_98_c62_n256))
