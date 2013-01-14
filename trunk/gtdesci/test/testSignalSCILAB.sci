function y = testSignalSCILAB(x)
    // Choose the signal parameters
    F = 8000;
    L = 200;
    // Evaluate the signal at the specified x's
    %v0_2 = size(x);    y = zeros(%v0_2(1),%v0_2(2));
    y = mtlb_i(y,mtlb_logic(x,">=",0),exp(-abs(L*mtlb_e(x,mtlb_logic(x,">=",0)))) .*sin(F*mtlb_e(x,mtlb_logic(x,">=",0))));
endfunction