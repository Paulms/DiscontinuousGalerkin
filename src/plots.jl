@recipe function f(sol::DGSolution; tidx = size(sol.t,1), uvars=0)
    xguide --> "x"
    yguide --> "u"
    labels = String[]
    for i in 1:size(sol.u[tidx],2)
      push!(labels,"u$i")
    end
    yvector = sol.u[tidx]
    ysvector = sol.uₛ[tidx]
    if uvars != 0
      yvector = sol.u[tidx][:,uvars]
      ysvector = sol.uₛ[tidx][:,uvars]
      labels = labels[uvars]
    end
    if typeof(labels) <: String
      label --> labels
    else
      label --> reshape(labels,1,length(labels))
    end
    @series begin
      seriestype  :=  :path
      sol.nodes, yvector
    end
    @series begin
      seriestype  :=  :scatter
      sol.face_nodes, ysvector
    end

end
