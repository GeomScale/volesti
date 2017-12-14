function [xopt,fval,lambda,exitflag,how]=mpt_solveLPi(f,A,B,Aeq,Beq,x0,lpsolver)

[xopt,fval,lambda,exitflag,how]=mpt_solveLP(f,A,B,Aeq,Beq);
