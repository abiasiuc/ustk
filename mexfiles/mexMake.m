function mexMake(file)
%MEXMAKE Compile MEX files
%
% Input
%	file		MEX file to compile
%				('gmmExtract', 'mstExtract', 'mstSmooth')

switch file
	case {'gmmExtract', 'gmmExtract.c'}
		mex mexfiles/gmmExtract.c ...
			mexfiles/utils/tqli.c ...
			mexfiles/utils/nrutil.c ...
			mexfiles/utils/tred2.c ...
			mexfiles/utils/matrix.c ...
			mexfiles/utils/gaussian.c ...
			mexfiles/utils/determineMicrostates.c ...
			-output mexfiles/mex_gmmExtract
		
	case {'mstExtract', 'mstExtract.c'}
		mex mexfiles/mstExtract.c ...
			mexfiles/utils/tqli.c ...
			mexfiles/utils/nrutil.c ...
			mexfiles/utils/tred2.c ...
			mexfiles/utils/matrix.c ...
			mexfiles/utils/gaussian.c ...
			mexfiles/utils/determineMicrostates.c ...
			-output mexfiles/mex_mstExtract
		
	case {'mstSmooth', 'mstSmooth.c'}
		mex mexfiles/mstSmooth.c ...
			mexfiles/utils/tqli.c ...
			mexfiles/utils/nrutil.c ...
			mexfiles/utils/tred2.c ...
			mexfiles/utils/matrix.c ...
			mexfiles/utils/gaussian.c ...
			mexfiles/utils/determineMicrostates.c ...
			-output mexfiles/mex_mstSmooth
		
	case {'gmmResponsibilities', 'gmmResponsibilities.c'}
		mex mexfiles/gmmResponsibilities.c ...
			mexfiles/utils/tqli.c ...
			mexfiles/utils/nrutil.c ...
			mexfiles/utils/tred2.c ...
			mexfiles/utils/matrix.c ...
			mexfiles/utils/gaussian.c ...
			mexfiles/utils/determineMicrostates.c ...
			-output mexfiles/mex_gmmResponsibilities
	
	case {'mstGaussian', 'mstGaussian.c'}
		mex mexfiles/mstGaussian.c ...
			mexfiles/utils/tqli.c ...
			mexfiles/utils/nrutil.c ...
			mexfiles/utils/tred2.c ...
			mexfiles/utils/matrix.c ...
			mexfiles/utils/gaussian.c ...
			-output mexfiles/mex_mstGaussian
	
	otherwise
		disp('Unknown mex file');
end

end

