% CPP_TEST_SCRIPT   Generate static library cpp_test from cpp_test.
% 
% Script generated from project 'cpp_test.prj' on 31-May-2024.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.CodeConfig'.
cfg = coder.config('lib','ecoder',false);
cfg.TargetLang = 'C++';
cfg.GenerateReport = true;
cfg.MaxIdLength = 1024;
cfg.ReportPotentialDifferences = false;
cfg.Toolchain = 'GNU gcc/g++ | CMake/gmake (64-bit Linux)';
cfg.GenCodeOnly = true;
cfg.BuildConfiguration = 'Release';

%% Define argument types for entry-point 'cpp_test'.
ARGS = cell(1,1);
ARGS{1} = cell(2,1);
ARGS{1}{1} = coder.typeof(0);
ARGS{1}{2} = coder.typeof(1i);

%% Invoke MATLAB Coder.
codegen -config cfg cpp_test -args ARGS{1}

