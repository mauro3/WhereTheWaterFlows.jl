# Original 37×48 grounding-line fixture, expressed as inactive runs per row.
function gl_order_fixture()
    mask = trues(37, 48)
    mask[6, 48:48] .= false
    mask[7, 47:48] .= false
    mask[8, 45:48] .= false
    mask[9, 45:48] .= false
    mask[10, 45:48] .= false
    mask[11, 45:48] .= false
    mask[12, 46:48] .= false
    mask[13, 43:48] .= false
    mask[14, 43:48] .= false
    mask[15, 42:48] .= false
    mask[16, 39:48] .= false
    mask[17, 39:48] .= false
    mask[18, 40:48] .= false
    mask[19, 43:48] .= false
    mask[20, [20:22; 35:36; 41:48]] .= false
    mask[21, [18:25; 35:40; 43:48]] .= false
    mask[22, [15:27; 29; 32:41; 44:48]] .= false
    mask[23, 13:48] .= false
    mask[24, 11:48] .= false
    mask[25, 11:48] .= false
    mask[26, 10:48] .= false
    mask[27, 10:48] .= false
    mask[28, 10:48] .= false
    mask[29, 9:48] .= false
    mask[30, 9:48] .= false
    mask[31, 5:48] .= false
    mask[32, 5:48] .= false
    mask[33, 4:48] .= false
    mask[34, 5:48] .= false
    mask[35, 5:48] .= false
    mask[36, 6:48] .= false
    mask[37, 6:48] .= false
    return (; mask, gl=findall(WWFS.get_boundary_cells(mask, true)))
end
