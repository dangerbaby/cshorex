
in.cd = 1;
in.Hrms = 1;
in.d = 2;
in.T = 10;
in.N = 10;
in.D = .01;
in.alpha = 0;
in.vegheight = 1;
in.U = undertow_linear(in.d,in.Hrms,in.T);
in.inonlin = 0;

[out0] = veg_stress_dissipation(in);
%[out1] = shear_stress(cf,Hrms,h,a,Tp,Ur,[0],1);
