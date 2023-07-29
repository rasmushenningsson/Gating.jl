mutable struct DataView
	data::DataMatrix
	reduced::Union{DataMatrix,Nothing}
	umapped::Union{DataMatrix,Nothing}
end
DataView(data) = DataView(data,nothing,nothing)


function get_reduced(v::DataView)
	if v.reduced === nothing
		@info "Computing SVD"
		v.reduced = svd(v.data; nsv=50)
	end
	v.reduced
end
function get_umapped(v::DataView)
	if v.umapped === nothing
		@info "Computing UMAP"
		v.umapped = umap(get_reduced(v), 2)
	end
	v.umapped
end


struct SingleCellGater <: AbstractGater
	data::DataMatrix # original data

	x::Observable{Any}
	y::Observable{Any}

	coords::Observable{Matrix{Float64}}
	# colors::Observable # TODO
	pl::Any

	poly_coords::Observable{Matrix{Float64}}
	poly_color::Observable{RGB{Float64}}
	poly_state::Base.RefValue{Symbol} # new, creating, finished
	poly_pl::Any

	stack::Vector{DataView}

	selection_x1::Base.RefValue{Union{Float64,Nothing}}
	selection_y1::Base.RefValue{Union{Float64,Nothing}}
end
function SingleCellGater(data::DataMatrix)
	coords = Observable(zeros(2,size(data,2)))

	pl = scatter(coords)
	x = pl.axis.xlabel
	y = pl.axis.ylabel

	poly_coords = Observable(zeros(2,3))
	poly_color = Observable{RGB{Float64}}(colorant"yellow")
	poly_pl = poly!(pl.axis, poly_coords; visible=true, strokewidth=2, color=poly_color, alpha=0.2)

	gater = SingleCellGater(data, x, y,
	                        coords,
	                        pl,
	                        poly_coords,
	                        poly_color,
	                        Ref(:finished),
	                        poly_pl,
	                        [DataView(data)],
	                        Ref{Union{Float64,Nothing}}(nothing),
	                        Ref{Union{Float64,Nothing}}(nothing),
	                       )

	on(x) do val
		x_changed(gater, val)
	end
	on(y) do val
		y_changed(gater, val)
	end
	on(events(pl.figure.scene).mousebutton; priority=10) do event
		mouse_handler(gater, event)
	end
	on(events(pl.figure.scene).keyboardbutton) do event
		key_handler(gater, event)
	end
	
	gater
end

curr_view(gater::SingleCellGater) = gater.stack[end]

function get_var(gater::SingleCellGater, name::String)
	v = curr_view(gater)

	if startswith(name, "pc")
		i = parse(Int, name[3:end])
		reduced = get_reduced(v)
		return reduced.matrix.Vt[i,:]
	elseif startswith(name, "umap")
		i = parse(Int, name[5:end])
		umapped = get_umapped(v)
		return umapped.matrix[i,:]
	else
		i = findfirst(isequal(name), v.data.var.name)
		ei = zeros(1,size(v.data,1))
		ei[i] = 1.0
		return vec(ei*v.data.matrix)
	end
end

function x_changed(gater::SingleCellGater, val)
	gater.coords[][1,:] .= get_var(gater, val)
	notify(gater.coords)
	# TODO: also update x axis name
end
function y_changed(gater::SingleCellGater, val)
	gater.coords[][2,:] .= get_var(gater, val)
	notify(gater.coords)
	# TODO: also update y axis name
end


function Base.push!(gater::SingleCellGater, ids::DataFrame)
	@assert gater.data.obs_id_cols == names(ids)

	ids = copy(ids; copycols=false)
	ids.__present__ .= true 

	# mask = gater.data.obs 
	merged = leftjoin(gater.data.obs, ids; on=gater.data.obs_id_cols, order=:left)
	mask = .!ismissing.(merged.__present__) 

	filtered = gater.data[:,mask]
	push!(gater.stack, DataView(filtered))
	gater.coords[] = vcat(get_var(gater, gater.x[])', get_var(gater, gater.y[])')

	gater
end

function Base.pop!(gater::SingleCellGater)
	@assert length(gater.stack)>1
	popped = pop!(gater.stack)

	gater.coords[] = vcat(get_var(gater, gater.x[])', get_var(gater, gater.y[])')

	select(popped.data.obs, popped.data.obs_id_cols)
end

function key_handler(gater::SingleCellGater, event)
	if event.key == Makie.Keyboard.left_shift
		if event.action == Makie.Keyboard.press
			gater.poly_state[] = :new
		else#if event.action == Makie.Keyboard.release
			gater.poly_state[] = :finished
			gater.poly_color[] = colorant"green"
		end
	end
end

function mouse_handler(gater::SingleCellGater, event)
	pl = gater.pl

	poly_state = gater.poly_state[]
	if poly_state != :finished
		if event.button == Mouse.left && event.action == Mouse.press
			poly_coords = poly_state==:new ? zeros(2,0) : gater.poly_coords[]
			mouse_coords = Float64[mouseposition(pl.axis.scene)...]
			poly_coords = hcat(poly_coords, mouse_coords)

			gater.poly_coords[] = poly_coords
			gater.poly_state[] = :creating
			gater.poly_color[] = colorant"yellow"

			return Consume()
		end 
	end

	# if Keyboard.left_shift in events(pl.figure.scene).keyboardstate
	# 	if event.button == Mouse.left && event.action == Mouse.press
	# 		gater.selection_x1[],gater.selection_y1[] = mouseposition(pl.axis.scene)
	# 		@info "selection started at ($(gater.selection_x1[]),$(gater.selection_y1[]))"

	# 		return Consume()
	# 	elseif event.button == Mouse.left && event.action == Mouse.release
	# 		x1,y1 = gater.selection_x1[], gater.selection_y1[] 
	# 		x2,y2 = mouseposition(pl.axis.scene)

	# 		x1,x2 = minmax(x1,x2)
	# 		y1,y2 = minmax(y1,y2)
	# 		@info "selection performed ($x1,$y1,$x2,$y2)"

	# 		mask = (x1 .<= gater.coords[][1,:] .<= x2) .& (y1 .<= gater.coords[][2,:] .<= y2)

	# 		if count(mask) > 0
	# 			@info "$(count(mask)) cells selected"
	# 			v = curr_view(gater)
	# 			ids = select(v.data.obs, v.data.obs_id_cols)
	# 			push!(gater, ids[mask,:])
	# 		else
	# 			@info "No cells selected, aborting"
	# 		end

	# 		return Consume()
	# 	end
	# end
end


function gate(data::DataMatrix)
	gater = SingleCellGater(data)
	display(gater.pl)
	gater
end


