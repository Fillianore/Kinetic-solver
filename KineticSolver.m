function KineticSolver()
% A kinetic solver code based on Matlab written by JianFu Chen
% under the guidance of PeiJun Hu and HaiFeng Wang on 2014/4/4
% at East China University of Science and Technology

% The reaction mechanism expressions:
% H2O(c) + #1 <-> OH_minus#1 + proton(c)
% OH_minus#1 + hole(c) <-> OH_rad#1
% OH_rad#1 <-> O_minus#1 + proton(c)
% O_minus#1 + O_minus#1 <-> O2_2minus#1 + #1
% O2_2minus#1 + hole(c) <-> O2_minus#1
% O2_minus#1 + hole(c) <-> O2(p) + #1
% Obr_2minus#2 + hole(c) <-> Obr_minus#2
% Obr_minus#2 + OH_rad#1 <-> ObrOH_minus#2 + #1
% ObrOH_minus#2 + hole(c) <-> ObrOH#2
% ObrOH#2 <-> ObrO_minus#2 + proton(c)
% ObrO_minus#2 + hole(c) <-> O2(p) + #2
% H2O(c) + #2 <-> H2O#2
% H2O#2 <-> OH_minus#2 + proton(c)
% OH_minus#2 <-> Obr_2minus#2 + proton(c)
%
ReactionNum = 14; % number of reactions
Q0 = [1 1];       % total free sites
T  = 300;         % reaction temperature : K
DigitNum = 100;   % digit number used in Mupad for modified Newton method
% Reaction reactant/production species concentration symbol
ReactionVarList = 'x_H2O C_hole x_O2 C_proton';
ReactionVarValue = [1     1     1E-7   1E-7   ];
% Reaction intermediate species coverage symbol, Q1_/Q2_ for site #1/#2
ReactionCovList1 = 'Q1_O2_2minus Q1_O2_minus Q1_OH_minus Q1_OH_rad Q1_O_minus Q1_v';
ReactionCovList2 = 'Q2_H2O Q2_OH_minus Q2_ObrOH Q2_ObrOH_minus Q2_ObrO_minus Q2_Obr_2minus Q2_Obr_minus Q2_v';
% The correspoding reaction free barrier energy (Ea) and free eneries (G0)
Ea = [ 0.5,     0,  0.35,  0.21,     0,     0,     0,  0.23,     0,   0.2,     0,     0,  0.34,  0.29];
G0 = [0.07, -0.17, -0.52, -1.31, -1.55, -1.13,  0.08, -1.59, -0.87, -0.14, -1.27, -0.96, -0.09, -0.28];
% set sample variable and value
SampleVar1 = 10.^(-10:0.5:0);% for reaction concentration sample : C_hole
ReactionVarSamplePosition = 2;
SampleVar2 = 0:0.1:1.2;    % for reaction barrier sample : Water dissociation/O-O coupling/hole diffusion
BarrierSamplePosition = [2,5,6,7,9,11]; % for hole diffusion, 1 for Water dissociation,[4 8] for O-O coupling
Path1Index = 6;
Path2Index = 11;
ResultFileName = 'Result.mat'; % file name of saved result
% set all input numeric values to be symbolic for high-precision numerical calculation
T = sym(T);
Q0 = sym(Q0);
Ea = sym(Ea);
G0 = sym(G0);
ReactionVarValue = sym(ReactionVarValue);
% ksymbol = [kf ; kr]; rate constant symbol matrix
kbasic = 'k';
for ir = 1:ReactionNum
    ixs = num2str(ir);
    order = ['(' ixs ')'];
    kf{ir} = [kbasic 'f' order]; %#ok
    kr{ir} = [kbasic 'r' order]; %#ok
end
kf = sym(mupadmex(kf));
kr = sym(mupadmex(kr));
ksymbol = [kf ; kr];
ReactionCovList = [ReactionCovList1 ' ' ReactionCovList2];
eval(['syms ' ReactionVarList ' ' ReactionCovList]);
% define the reaction forward and reverse rate equations
ForwardRateEquation = [ kf(1)*x_H2O*Q1_v;...
    kf(2)*C_hole*Q1_OH_minus;...
    kf(3)*Q1_OH_rad;...
    kf(4)*Q1_O_minus^2;...
    kf(5)*C_hole*Q1_O2_2minus;...
    kf(6)*C_hole*Q1_O2_minus;...
    kf(7)*C_hole*Q2_Obr_2minus;...
    kf(8)*Q1_OH_rad*Q2_Obr_minus;...
    kf(9)*C_hole*Q2_ObrOH_minus;...
    kf(10)*Q2_ObrOH;...
    kf(11)*C_hole*Q2_ObrO_minus;...
    kf(12)*x_H2O*Q2_v;...
    kf(13)*Q2_H2O;...
    kf(14)*Q2_OH_minus];
ReverseRateEquation = [ kr(1)*C_proton*Q1_OH_minus;...
    kr(2)*Q1_OH_rad;...
    kr(3)*C_proton*Q1_O_minus;...
    kr(4)*Q1_v*Q1_O2_2minus;...
    kr(5)*Q1_O2_minus;...
    kr(6)*x_O2*Q1_v;...
    kr(7)*Q2_Obr_minus;...
    kr(8)*Q1_v*Q2_ObrOH_minus;...
    kr(9)*Q2_ObrOH;...
    kr(10)*C_proton*Q2_ObrO_minus;...
    kr(11)*x_O2*Q2_v;...
    kr(12)*Q2_H2O;...
    kr(13)*C_proton*Q2_OH_minus;...
    kr(14)*C_proton*Q2_Obr_2minus];
RateEquation = ForwardRateEquation - ReverseRateEquation;
ReversibilityEquation = ReverseRateEquation./ForwardRateEquation;
% define the reaction rate differential equations
% dQ1_O2_2minus/dt = r(4) - r(5)
% dQ1_O2_minus/dt = r(5) - r(6)
% dQ1_OH_minus/dt = r(1) - r(2)
% dQ1_OH_rad/dt = r(2) - r(3) - r(8)
% dQ1_O_minus/dt = r(3) - 2*r(4)
% dQ1_v/dt = r(4) - r(1) + r(6) + r(8)
% dQ2_H2O/dt = r(12) - r(13)
% dQ2_OH_minus/dt = r(13) - r(14)
% dQ2_ObrOH/dt = r(9) - r(10)
% dQ2_ObrOH_minus/dt = r(8) - r(9)
% dQ2_ObrO_minus/dt = r(10) - r(11)
% dQ2_Obr_2minus/dt = r(14) - r(7)
% dQ2_Obr_minus/dt = r(7) - r(8)
% dQ2_v/dt = r(11) - r(12)
SteadyStateEquation(1)  = RateEquation(4) - RateEquation(5);
SteadyStateEquation(2)  = RateEquation(5) - RateEquation(6);
SteadyStateEquation(3)  = RateEquation(1) - RateEquation(2);
SteadyStateEquation(4)  = RateEquation(2) - RateEquation(3) - RateEquation(8);
SteadyStateEquation(5)  = RateEquation(3) - 2*RateEquation(4);
SteadyStateEquation(6)  = RateEquation(4) - RateEquation(1) + RateEquation(6) + RateEquation(8);
SteadyStateEquation(7)  = RateEquation(12) - RateEquation(13);
SteadyStateEquation(8)  = RateEquation(13) - RateEquation(14);
SteadyStateEquation(9)  = RateEquation(9) - RateEquation(10);
SteadyStateEquation(10) = RateEquation(8) - RateEquation(9);
SteadyStateEquation(11) = RateEquation(10) - RateEquation(11);
SteadyStateEquation(12) = RateEquation(14) - RateEquation(7);
SteadyStateEquation(13) = RateEquation(7) - RateEquation(8);
SteadyStateEquation(14) = RateEquation(11) - RateEquation(12);

ReactionVarString = regexp(ReactionVarList,' ','split');
ReactionVarsymbol = sym(mupadmex(ReactionVarString));
ReactionCovString1 = regexp(ReactionCovList1,' ','split');
ReactionCovsymbol1 = sym(mupadmex(ReactionCovString1));
ReactionCovString2 = regexp(ReactionCovList2,' ','split');
ReactionCovsymbol2 = sym(mupadmex(ReactionCovString2));
ReactionCovsymbol = [ReactionCovsymbol1 ReactionCovsymbol2];
CoverageNum = length(ReactionCovsymbol);
% The mass balances of free sites
SteadyStateEquation(6)  = sum(ReactionCovsymbol1) - Q0(1);
SteadyStateEquation(14) = sum(ReactionCovsymbol2) - Q0(2);

kB = sym('1.3806505e-23');  % Boltzmann constant
h = sym('6.62606957e-34');  % Plank constant
e = sym('1.60217653e-19');  % elementary charge
SampleNum1 = length(SampleVar1);
SampleNum2 = length(SampleVar2);
SampleNum = SampleNum1*SampleNum2;
ReactionRate = zeros(SampleNum,ReactionNum);
Reversibility = zeros(SampleNum,ReactionNum);
CoverageSolution = zeros(SampleNum,CoverageNum);
% solve the kinetic
for ix = 1:SampleNum1
    ReactionVarValue(ReactionVarSamplePosition) = SampleVar1(ix);
    for jx = 1:SampleNum2
        Ea(BarrierSamplePosition) = SampleVar2(jx);
        index = (ix - 1)*SampleNum2 + jx;
        kvalue = RateConstantCalculation(Ea,G0,T,kB,e,h);
        [ReactionRate(index,:),Reversibility(index,:),CoverageSolution(index,:)] = NumericSolve(index,RateEquation...
            ,ReversibilityEquation,SteadyStateEquation,ReactionVarsymbol,...
            ReactionVarValue,ksymbol,kvalue,ReactionCovsymbol,DigitNum);
    end
end
% save the calculation result
save(ResultFileName,'ReactionRate','Reversibility','CoverageSolution');
% plot the Turnover frequency
TurnoverFreq1 = reshape(ReactionRate(:,Path1Index),SampleNum2,SampleNum1);
TurnoverFreq2 = reshape(ReactionRate(:,Path2Index),SampleNum2,SampleNum1);
TurnoverFreq = TurnoverFreq1 + TurnoverFreq2;
% Selectivity = TurnoverFreq1./TurnoverFreq2;
[XSample,YSample] = meshgrid(SampleVar1,SampleVar2);
contourf(XSample,YSample,log10(TurnoverFreq),100,'LineStyle','none');colorbar;colormap Jet
xlabel('log_{10}(C_{hole})');set(gca,'xscale','log'); caxis([-19 4.5]);
ylabel('Ea_{hole}(eV)');
title('log_{10}(Turnover frequency)');snapnow;
% surf(XSample,YSample,Selectivity,'LineStyle','none');colorbar;colormap Jet
% view(2)
% xlabel('log_{10}(C_{hole})');set(gca,'xscale','log');
% ylabel('Ea_{hole}(eV)');
% title('log_{10}(Selectivity)')
function [ReactionRate,Reversibility,CoverageSolution] = NumericSolve...
    (index,RateEquation,ReversibilityEquation,SteadyStateEquation,ReactionVarsymbol...
    ,ReactionVarValue,ksymbol,kvalue,ReactionCovsymbol,DigitNum)
% eval the SteadyStateEquation
SteadyStateEquation1 = subs(subs(SteadyStateEquation,ReactionVarsymbol,ReactionVarValue));
SteadyStateEquation2 = subs(SteadyStateEquation1,ksymbol,kvalue);
% transform SteadyStateEquation to strings
SteadyStateEquationStrings = char(SteadyStateEquation2);
SteadyStateEquationStrings = SteadyStateEquationStrings(10:end - 3);
ReactionCovStrings = char(ReactionCovsymbol);
ReactionCovStrings = ReactionCovStrings(10:end - 3);
% send the SteadyStateEquation to mupad, and solved by modified Newton method
ReactionSolution = evalin(symengine,['DIGITS := ' num2str(DigitNum) ': numeric::fsolve({', SteadyStateEquationStrings '},[' ReactionCovStrings '])']);
CoverageSolution = mygetsol(ReactionSolution,ReactionCovStrings);
fprintf(['\nRun the Sample ' num2str(index) ' now\n']);
disp('Coverages:');
for ir = 1:length(ReactionCovsymbol)
    disp([ char(ReactionCovsymbol(ir)) ' = ' num2str(double(CoverageSolution(ir)),'%16.10e')]);
end
ReactionRate = subs(RateEquation,ReactionVarsymbol,ReactionVarValue);
ReactionRate = subs(subs(ReactionRate,ksymbol,kvalue),ReactionCovsymbol,CoverageSolution);
Reversibility = subs(ReversibilityEquation,ReactionVarsymbol,ReactionVarValue);
Reversibility = subs(subs(Reversibility,ksymbol,kvalue),ReactionCovsymbol,CoverageSolution);
fprintf('\n            Rate                  Reversibility\n');
ReactionNum = length(ReactionRate);
ReactionIndex = sym('R%d',[ReactionNum,1]);
for ir = 1:ReactionNum
    fprintf('%-8s:   %s\n', char(ReactionIndex(ir)), num2str(double([ReactionRate(ir) Reversibility(ir)]),'%22.10e'));
end
function kvalue = RateConstantCalculation(Ea,G0,T,kB,e,h)
% calculate the forward and reverse reaction rates
% based on free barrier energy (Ea) and free eneries (G0)
keV = e/T/kB;
kBT_h = kB/h*T;
Keq = exp(-G0*keV);
kf = exp(-Ea*keV)*kBT_h;
kr = kf./Keq;
kvalue = [kf ; kr];
function CoverageSolution = mygetsol(Reactionsolution,ReactionCovsymbol)
% get the solution from the symbol matrix
solutions = strrep(strrep(char(Reactionsolution),'matrix([',''),']','');
eval([strrep(strrep(solutions(2:end - 1),' == ','_new = sym('''),',',''');') ''');']);
CoverageSolution = eval(['['  [strrep(ReactionCovsymbol,',','_new, ') '_new'] ']']);
