function mpt_demo_PowerNetworkSystem
%
% Formulation of MPC problem for power network system
%

close all

% load data provided by S. Riverso <stefano.riverso@unipv.it>, University of Pavia
load dataSim

% construct discrete-time LTI model
% ( The prediction model contains measurable disturbance deltaPload which cannot
% be handled directly with the current modeling approach. Alternatively, we
% consider augmented model where the measurable disturbance is a part of
% dynamical matrix and formulate MPC problem not to penalize these additional states.)  
model = LTISystem('A',[A,L;zeros(4,16),eye(4)],'B',[B; zeros(4)],'C',[C,zeros(16,4)],'D',D,'Ts',Ts);

% simulation model x+ = A*x + B*u + L*deltaPload
simulation_model = LTISystem('A',A,'B',[B, L],'C',C,'D',[D, zeros(16,4)],'Ts',Ts);

%% formulate MPC problem 
ctrl = MPCController(model);

% set the prediction horizon
ctrl.N = 15;

% input constraints
ctrl.model.u.min = [-0.5; -0.65; -0.65; -0.55];
ctrl.model.u.max = [0.5; 0.65; 0.65; 0.55];

% state constraints
ctrl.model.x.min = [-0.1; -Inf(19,1)];
ctrl.model.x.max = [0.1; Inf(19,1)];
 
% add reference signals
ctrl.model.x.with('reference');
ctrl.model.u.with('reference');

% zero terminal set
ctrl.model.x.with('terminalSet');
% ctrl.model.x.with('terminalPenalty');
% ctrl.model.x.terminalPenalty = QuadFunction(eye(16));

% add quadratic penalty on the outputs
ctrl.model.u.penalty = QuadFunction(R);
% penalty must be positive definite
ctrl.model.x.penalty = QuadFunction(blkdiag(Q,zeros(4)));


%% Simulate the closed loop
x0 = zeros(16,1);
X = [];
U = [];
Y = [];
N_sim = length(deltaPref.signals.values);
simulation_model.initialize(x0);
eta = 0;
phi = 0;
x = zeros(16,1);
u = zeros(4,1);
for k=1:N_sim

    % time varying reference for inputs
    ctrl.model.u.reference = uO(k,:)';

    % time varying reference for states
    dPl = deltaPload.signals.values(k,:)';
    ctrl.model.x.reference = [xO(k,:)'; dPl];
    ctrl.model.x.terminalSet = Polyhedron('Ae',eye(20),'be',ctrl.model.x.reference);    

    % evaluate the performace criteria before state update
    eta = eta + (x(1:16)-ctrl.model.x.reference(1:16))'*Q*(x(1:16)-ctrl.model.x.reference(1:16)) + (u-dPl)'*R*(u-dPl);
    Ptie12 = 4*(x(1)-x(5));
    Ptie23 = 2*(x(5)-x(9));
    Ptie34 = 2*(x(9)-x(13));
    phi = phi + sum(abs(Ptie12)*Ts + abs(Ptie23)*Ts + abs(Ptie34)*Ts);
    
    % evaluate MPC
    u = ctrl.evaluate([x;dPl]);
    
    % load precomputed inputs from stored data
    %u = deltaPref.signals.values(k,:)';
    
    % update the process
    [x,y] = simulation_model.update([u;dPl]);
    
    % store the signals
    U = [U, u(1:4)];
    X = [X, x(1:16)];
    Y = [Y, y];
    fprintf('Simulating... %d/%d\n',k,N_sim);

end

% evaluate the performance criteria
eta = eta/N_sim
phi = phi/N_sim

% plot the state trajectories
figure
plot(1:N_sim,X,'LineWidth',2)
xlabel('Simulation steps')
title('States')

% plot the control inputs
figure
plot(1:N_sim,U,'LineWidth',2)
xlabel('Simulation steps')
title('Control inputs')

end
