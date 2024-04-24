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
  coords = coordinates.(verts)
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


