%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script for counting  arithmetic operations in various compact models to
% compare complexity!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MOD  = BSIM3v3_2_4_ModSpec();

% Find out the number of unknowns to be supplied to the model
n_vx = length(MOD.OtherIONames(MOD));
vx   = rand(n_vx, 1);

n_vy = length(MOD.InternalUnkNames(MOD));
vy   = rand(n_vy, 1);

n_vu = length(MOD.uNames(MOD));
vu   = rand(n_vu, 1);

% Construct a `flag` to be passed into the model
flag    = struct();
flag.fe = 1;
flag.qe = 1;
flag.fi = 1;
flag.qi = 1;

% Dummy line to enable/disable profining
set_debug = 1;

% Evaluate the core function
[fe, qe, fi, qi] = MOD.fqei(vx, vy, vu, flag, MOD);
