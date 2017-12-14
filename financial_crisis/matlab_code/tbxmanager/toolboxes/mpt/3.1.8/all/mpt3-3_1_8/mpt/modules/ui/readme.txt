User Interface for Modeling and Control Synthesis for MPT 3.0
=============================================================

The proposed user interface (UI) aims to fix the design flaws of the
old MPT 2.6.x interface. The main points behind the new UI are:
* consistency with respect to working with various prediction models
(LTI, PWA, MLD, compositions of various subsystems)
* consistency of mathematical formulation between on-line and explicit
controllers
* easy extensibility without the need to maintain a single monolithic
implementation
* ability to seamlessly import MPT 2.6.x formulations
* promotion of the "correct" workflow where tuning is done on on-line
controllers first, only converting them to the explicit form when the
tuning is complete (avoiding frequent re-computation of the explicit
solution) 

The UI consists of three primary layers, described in more details below.
1) Classes representing prediction models
2) Classes representing controllers
3) Class representing a generic closed-loop system

1. Representation of prediction models
--------------------------------------

The UI allows to define four basic types of prediction models:
1) Linear Time Invariant (LTI) models
2) Piecewise Affine (PWA) models
3) Mixed-Logical Dynamical (MLD) models
4) Compositions consisting of several interconnected LTI/PWA/MLD
models

1.1 LTI models
--------------

LTI models describe linear (in fact, affine) time-invariant systems
of the form
  x^+ = A x + B u + f
    y = C x + D u + g

Such models, represented by the @LTISystem class, can be constructed
in three ways:

1) By providing the matrices:
      lti = LTISystem('A', A, 'B', B, 'C', C, 'D', D, 'f', f, 'g', g)
   The 'f' and 'g' parameters can be omitted.

2) By importing from state-space models of the Control Toolbox:
      ss = ss(A, B, C, D, Ts)
      lti = LTISystem(ss)
   Note that only discrete-time state-space objects can be imported.

3) By importing from MPT2 sysStruct:
      lti = LTISystem(sysStruct)

The returned object contains:
1) matrices A, B, C, D, f, g of the state-update and output equations
   (read-only) 
2) dimensions nx, nu, ny (read-only)
3) sampling time Ts (read-only)
4) signals representing states (lti.x), inputs (lti.u) and outputs
   (lti.y) of the system (read/write)
5) domain of the system in the state-input space (can be set using the
   setDomain method)

Properties of these signals (e.g. constraints and penalties) can be
defined as described in Section 1.4.

The class implements following methods:
* initialize(x0): Sets the current state of the system
* getStates(): Returns the current state of the system
* update(u): Performs the state update using the control action "u"
* output(u): Computes the output (note that providing the control
  action is only needed if the system contains direct feed-through)
* setDomain(type, D): sets the domain of the system either in the state 
  space (type='x'), input space (type='u') or state-input space (type='xu')
* plot(): to visualize the optimized open-loop profiles of states,
  inputs and outputs (requires prior call to an optimization routine,
  explained in Section 2.x)

The object remembers the value of its state and updates it whenever
update() is called. This means that evolution of the system can be
obtained in the following way:
    lti.initialize(x0)
    lti.update(u0)
    lti.update(u1)
    lti.update(u2)
    lti.getStates()

See mpt3_lti_demo1.m, mpt3_lti_demo2.m, mpt3_lti_demo3.m, 
mpt3_lti_demo4.m, mpt3_lti_demo5.m for examples.

1.2 PWA models
--------------

PWA models encode the PWA dynamics

x^+ = A_i*x + B_i*u + f_i  IF  [x; u] \in Dom_i
  y = C_i*x + D_i*u + g_i  IF  [x; u] \in Dom_i

PWA models can be created in three ways:

1) By providing an array of LTI systems, each of which represents one mode
   of the PWA system:

      dyn1 = LTISystem(...)
      dyn1.setDomain('x', Polyhedron(...))  % Region of validity

      dyn2 = LTISystem(...)
      dyn2.setDomain('x', Polyhedron(...))  % Region of validity

      pwa = PWASystem([dyn1, dyn2])

2) By converting an MLD description to a PWA form:

      mld = MLDSystem(...)
      pwa = PWASystem(mld)

3) By importing from MPT2's sysStruct:

      opt_sincos
      pwa = PWASystem(sysStruct)

Internally, PWA models are represented by the @PWASystem class, which
contains following fields (all read-only):
* A, B, C, D, f, g: cell arrays of matrices
* domain: array of polyhedra denoting region of validity of corresponding
  local affine models
* nx, nu, ny, ndyn: number of states, inputs, outputs and dynamics
* Ts: sampling time
* x, u, y, d: signals representing states, inputs, outputs and binary
  mode selectors

Properties of signals 'x', 'u', 'y', 'd' (e.g. constraints and penalties) 
can be defined as described in Section 1.5.

The class implements the same methods as in the LTI case (Section 1.1).

1.3 MLD models
--------------

MLD models are represented by the @MLDSystem class and can be imported
from a HYSDEL source file by calling

     mld = MLDSystem('filename.hys')

If the source file contains symbolic parameters, their values can be
provided as a structure in the second input argument:

     mld = MLDSystem('filename.hys', struct('param1', value1))

MLD models can be converted to a PWA form using the toPWA() method:

     pwa = mld.toPWA();

The class implements the same methods as in the LTI case (Section 1.1).

1.4 Signals
-----------

As mentioned in Section 1.1-1.3, each prediction model contains
"signals", which represent the state, input and output variables of
the model. 

A "signal", represented by the @SystemSignal class, is a basic
primitive which does two things:
1) it represents the prediction of a given quantity over a prediction
   horizon as a YALMIP variable 
2) allows to define basic properties of the quantity

A "signal" is created by calling the s=SystemSignal(n)
constructor where "n" is the dimension of the signal. 
(NOTE: the user never has to create signals manually.)
Once the signal is created, following basic properties can be set by
the user:
* lower bound (s.min = [...])
* upper bound (s.max = [...])
* penalty function (s.penalty = ...)

As a specific example, consider the following LTI model of a double
integrator:

    Ts = 1
    s = c2d(ss(tf(1, [1 0 0])), Ts)
    lti = LTISystem(s)

Then the "lti" object contains following signals:

* lti.x: the 2x1 state vector
* lti.u: the 1x1 input
* lti.y: the 2x1 output vector

To set the lower/upper bound on the input, one would call

    lti.u.min = -1;
    lti.u.max = 1;
 
To set a quadratic state penalty of the form x'*Q*x on the state with
Q=eye(2), call

    lti.x.penalty = Penalty(eye(2), 2) 

NOTE: Definition of penalties is described in Section 1.5.

Besides the minimal set of properties described above, arbitrary
additional properties can be added on-the-fly by the concept of custom
properties, called "filters". Examples include polyhedral constraints, soft
constraints, or move blocking. These custom properties are described by 
separate methods of the SystemSignal class, prefixed by "filter_" (see 
"methods SystemSignal").

Examples:
* polyhedral set constraint (for states, inputs, outputs):

    lti.x.with('setConstraint') % enable the property
    lti.x.setConstraint = Polyhedron(...);  % specification of the property

* marking variables as binary (for states, inputs, outputs):

    lti.u.with('binary')
    lti.u.binary = true;

* move-blocking constraint u(1)=u(2)=u(3):

    lti.u.with('block')
    lti.u.block.from = 1
    lti.u.block.to = 3

* soft constraints (for states, inputs, outputs):

    lti.x.with('softMin')  % soft lower bound
    lti.x.with('softMax')  % soft upper bound

* penalty on the final predicted state:

    lti.x.with('terminalPenalty')
    lti.x.terminalPenalty = Penalty(...)

* polyhedral constraints on the final predicted state:

    lti.x.with('terminalSet')
    lti.x.terminalSet = Polyhedron(...)

* adding a reference trajectory for MPC (for states, inputs, outputs):

    lti.x.with('reference')
    lti.x.reference = [...]

* penalty on slew rate (for states, inputs, outputs):

    lti.u.with('deltaPenalty')
    lti.u.deltaPenalty = Penalty(...)

NOTE: only a single instance of each property can be added to the signal.

Since the framework of custom properties is general, it is perfectly
fine to apply blocking and delta-penalties to arbitrary variables, not
just to inputs. To minimize slew rate of outputs, one would set

    lti.y.with('deltaPenalty')
    lti.y.deltaPenalty = Penalty(...)

Keeping the custom properties as separate methods allows for great
flexibility as new properties can be added at any time.

1.5 Penalties
-------------

The @Penalty class defines penalty functions which are applied to
individual signals when forming the MPC cost function. Each signal can
be assigned different penalty function (and time-varying penalties can
be achieved by the concept of "custom properties" as discussed in
Section 1.5).

Penalty can be created by calling

    p = Penalty(Q, n)

which describes ||Q*z||_n, where "Q" is the penalty matrix, "z" is the 
variable to be penalized, and "n" is the desired norm:
    * n=0 leads to the penalty Q*z (ideal for maximization of physical
    quantities where "z" only takes non-negative values)
    * n=1, n=Inf give the standard 1/Inf norm
    * n=2 is the squared 2-norm z'*Q*z

To set a penalty on a signal, simply assign the penalty to the
signal's "penalty" property, e.g.:

    lti.x.penalty = Penalty(eye(2), 2)

Note that, unlike MPT 2.6.x, the new UI allows to use different norms
for different signals. I.e. it is perfectly fine to define

    lti.x.penalty = Penalty(eye(2), 2) % quadratic penalty on states
    lti.u.penalty = Penalty(1, Inf) % Inf-norm penalization of inputs
    lti.y.penalty = Penalty(-1, 0) % maximize the output

The penalty object can be evaluated for a given value of the argument
"z" by calling

    value = penalty.evaluate(z)

The penalization matrix "Q" and the norm can be changed on-the fly by 
directly modifying the "penalty.Q" and "penalty.norm" fields:

    lti.x.penalty = Penalty(eye(2), 2)
    lti.x.penalty.Q = 100*eye(2)  % amplify the penalization
    lti.x.penalty.norm = 1;       % change to L1-norm

When building the cost function, penalty terms are applied to the
signals at each step of the prediction horizon k=0...N-1. Note that no
penalty is applied on the final predicted state (k=N) for consistency
reasons. The terminal-state penalty can be added on-demand by applying
the "terminalPenalty" custom property:

    lti.x.with('terminalPenalty')
    lti.x.terminalPenalty = Penalty(...)

2. Representation of controllers
--------------------------------

The UI promotes the workflow where we always start with an on-line MPC
controller. This one is created upon calling the MPCController()
constructor as follows:
    controller = MPCController(model, prediction_horizon)
where "model" can be LTI/PWA/MLD/Composition. The returned object
contains two properties:
* model: represents the prediction model
*     N: prediction horizon (read/write)

The most sane workflow would look like this:
    model = LTISystem(...)
    controller = MPCController(model, horizon)
    controller.model.x.min = ...
    controller.model.x.max = ...
    controller.model.x.penalty = ...
and so on. 

When the controller is fully specified, it does just one thing --
computes the optimal control inputs for a specified value of the
initial state condition by calling

    u = controller.evaluate(x0)

By default, the evaluate() method returns the first element of the 
optimized sequence of control inputs. The open-loop optimizer can be
obtained by 

    U = controller.evaluate(x0, 'full')

If the problem is infeasible, the returned optimizer "U" will be a NaN
matrix and "J" will be Inf.

Once the optimization was performed, one can visualize the open-loop
profiles by calling

    controller.model.plot()

This will plot the profiles of optimal states, inputs and
outputs. Numerical values of the optimal profiles can be obtained by

    x_openloop = controller.model.x.value()
    u_openloop = controller.model.u.value()
    y_openloop = controller.model.y.value()

If the controller performs as desired, its explicit form can be
computed by calling

    explicit_controller = controller.toExplicit()

This creates an instance of the @EMPCController class, which behaves
in exactly the same fashion as the on-line controller.

The explicit controller supports the evaluate() method which returns
same outputs as in the on-line MPC case. 

Explicit controllers contain following read-only properties:
* partition: polyhedral partition of the optimizer
*  feedback: PolyUnion representation of the feedback law
*      cost: PolyUnion representation of the cost function

Moreover, the full optimizer (which combines the partition, feedback and
cost) is stored in the "optimizer" property.

Explicit controllers can be used to obtain values of control moves using
the evaluate() method:

    u = controller.evaluate(x0) % only the first control input
    U = controller.evaluate(x0, 'full') % open-loop optimizer

The polyhedral partition of an explicit controller can be plotted by
directly calling the plot() method on the partition:

    controller.partition.plot()

To plot the feedback law and the cost, call fplot() on the corresponding
properties:

    controller.feedback.fplot()
    controller.cost.fplot()

From the user's perspective, there is no difference whether he works
with on-line or with explicit controllers -- both use the same
mathematical setup and provide the same output (unlike MPT 2.6.x where
explicit controller in certain scenarios used different setups and
only provided closed-loop control action).

The following example illustrates how to perform a closed-loop
simulation over a fixed number of simulation steps with a different
simulation model:

    prediction_model = LTISystem(...)
    simulation_model = PWAModel(...)
    controller = MPCController(prediction_model, horizon)
    simulation_model.initialize(x0)
    for k = 1:N_steps
      u = controller.evaluate(x0)
      simulation_model.update(u)
      x0 = simulation_model.getStates()
    end

Organization of the control module allows it to easily support
other type of controllers, such as linear state-feedback or heuristic
controllers defined by PWA/MLD systems. All that needs to be done is
to derive from the @AbstractController and the new controller class
has to implement the evaluate() method.

3. Formal notion of closed-loop systems
---------------------------------------

The main drawback of MPT 2.6.x is that it has a flawed concept of
closed-loop systems, which are solely represented by the controller
object. While this simplifies closed-loop simulation, this approach
falls short in any other scenario, such as analysis (invariance,
stability). Moreover, in MPT 2.6.x it is extremely hard to define
non-MPC controllers, such as LQR state-feedback, and to analyze them.

Therefore the new UI allows to define a closed-loop system consisting
of a controller (does not matter which type), along with the
simulation model:
    loop = ClosedLoop(controller, model)

In the current version, the class implements just a single method to
perform closed-loop simulations:

    sim = loop.simulate(x0, N_steps)

It returns a structure with fields "X", "U", and "Y", containing the
closed-loop profiles of states, inputs and outputs, respectively.

Once the polyhedral library is finished, we will start implementing
the analysis functions on closed-loop systems, including invariance
analysis, reachable set computation and construction of Lyapunov
functions. Needless to say, all these tasks will require that the
controller is an explicit one.

See mpt3_lti_demo3.m
