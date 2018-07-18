function cart0 = Cart2AZ(cart)
    % Convert from cartesian to azimuth
    cart0 = reorient(cart)
    
    ndx1 = cart >= 0. & cart < 90.
    cart0(ndx1) = 90. - cart(ndx1)    
    ndx2 = cart >= 90. & cart < 180.
    cart0(ndx2) = 270. - cart(ndx2)
    ndx3 = cart >= 180. & cart < 360.
    cart0(ndx3) = 360. - cart(ndx3) + 90.
end

function arg = reorient(arg)
    ndx = arg < 0.0;
    arg(ndx) = arg(ndx) + 360.
end

