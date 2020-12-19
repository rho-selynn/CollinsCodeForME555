%//Condenser Sizing//

%input values below:
Nexch = 6;
ms = 55.56/Nexch; %steam mass flow rate/# of condensers
Rfi = 0.00015; Rfo = 0.00013; %fouling factors
Tsat = 41.51; Tmo = 39.51; Tmi = 15; %temperature values
hin = 1908.132; hout = 173.88; %enthalpy values entering and exiting condenser
u_cw = 2; %velocity of cooling water

%Cooling water properties @ 27.25C

k = 0.614; v = 8.506e-7; Prw = 5.778; p = 996.8; Cpw = 4.1789; u = 8.479e-4;

%Steam properties @ 8kPa and Tsat = 41.51

hf = 173.88; hfg = 2403.10; Prs = 4.191; pl = 991.214; Cps = 4.1789; ul = 635.4e-6; kl = 0.634; 

%properties of tubes (1-1/4in Schedule 40 square pitch):
 
Di = 0.035052; Di_inches = 1.38; Do = 0.04216; Do_inches = 1.660; tw = 0.003556; tw_inches = 0.140; PT = 0.03969;  
Np = 6; %Number of passes in the shell
kw = 48.9; %thermal conductivity of 1 Cr–V(0.2% C, 1.02% Cr, 0.15% V) alloy steel

%Start code:

Tlm = ((Tsat - Tmo) - (Tsat - Tmi)) / log((Tsat - Tmo) / (Tsat - Tmi));                                                                                 

Tm = (Tmi + Tmo)/2;

Q = ms*(hin - hout); %required rate of heat transfer from steam to cooling water in kW

mcw = Q / (Cpw*(Tmo -Tmi)); %mass flow rate of the cooling water in kg/s

Ntube = 4*mcw / (u_cw*p*pi*Di^2); %number of tubes/pass
Ntube = round(Ntube, 0);
Nr = Ntube*Np; %number of tubes required

N = 1008; Npass = 1008/6; Ds = 1.524; Ds_inches = 60; %shell tube layouts TEMA standard

ucw = 4*mcw / (p*pi*Di^2*Npass) %velocity of cooling water inside condenser tubes for the shell diameter
ReD = ucw*p*Di / u; %tube side Reynolds number
ReD = round(ReD, 0);

f = (0.790*log(ReD)- 1.64)^-2; %frictional factor for Nuesslet number equation below
NuD = ((f/8)*(ReD - 1000)*Prs) / (1 + ((12.7*(f/8)^0.5) * (Prs^(2/3) - 1))); %Nusselt number for a smooth pipe

hi = NuD*(k/Di); %tube side convection coefficient

NT = Ds/PT; %average number of transverse tubes
NT = round(NT, 0);

Dm = (Do - Di) / (log(Do/Di));  
Rt = Rfo + ((1+hi*Rfi)/hi)*(Do/Di) + (tw/kw)*(Do/Dm); %thermal resistance (m^2*K/W)

dTin = Tsat - Tmi;
dTout = Tsat - Tmo;
Tfoi = 5.291; Tfoo = 0.198; %guesses for fouled surface at the inlet and outlet of the tube (C or K)

hfgi = hfg + 0.68*Cpw*(Tfoi); %tube inlet (kJ/kg)
hfgo = hfg + 0.68*Cpw*(Tfoo); %tube outlet (kJ/kg)

hoi = 0.729*((pl^2*9.81*hfgi*10^3*kl^3)/(ul*Tfoi*Do))^0.25 * (1/NT^(1/6)); %shell side convection coefficient at inlet (W/m^2*K)
hoo = 0.729*((pl^2*9.81*hfgo*10^3*kl^3)/(ul*Tfoo*Do))^0.25 * (1/NT^(1/6)); %shell side convection coeffecient at outlet (W/m^2*K)

%Overall heat transfer Coefficients:

Uwi = (Rt + 1/hoi)^-1; %inlet (W/m^2*K)
Uwo = (Rt + 1/hoo)^-1; %outlet (W/m^2*K)

%Test initial guesses for Tfoi and Tfoo:

Tfoi = dTin*(1-Rt*Uwi);
Tfoo = dTout*(1-Rt*Uwo);

U = (Uwi + Uwo)/2; %overall convection coefficient is average of inlet and outlet

Ao = Q*10^3 / (U*Tlm); %required surface area in m^2

%Length of Condenser:

L = Ao / (Np*Npass*pi*Do); % length in meters
Lft = (L*39.3701)/12; %length in feet

%Tube side pressure drop:

f = 0.046*ReD^-0.2; %friction factor
G = p*ucw; %mass velocity
Pt = 4*f*((L*Np)/Di)*((G^2)/(2*p)); %pressure drop in the tubes(Pa)
Pr = 4*Np*((p*ucw^2)/2); %pressure drop in the return fitting (Pa)
Ptot = (Pt + Pr)/10^3; %total pressure drop through the tubes (kPa)
Ptot_psi = Ptot*0.145038; %convert to psi

%Pumping power required for the cooling water:

PWR = ((mcw*Ptot)/(p*0.85)); %power in kW

%Display sizing values:

fprintf('Final Design Summary: \n\n\n')
fprintf('Number of Heat Exchangers  %d \n', Nexch)
fprintf('Shell Diameter             %.3f m (max 60 in) \n', Ds)
fprintf('Shell Length               %.2f m (max 35ft) \n', L)
fprintf('Inner Diameter             %.4f m \n', Di)
fprintf('Outer Diameter             %.4f m \n', Do)
fprintf('Tube Pressure Drop         %.3f kpa \n', Ptot)
fprintf('Pumping Power              %.2f kW \n', PWR)
fprintf('Number of Tubes Required   %d tubes \n', Nr)
fprintf('Transverse Tubes           %d tubes \n', NT)
fprintf('Number of Tubes            %d tubes \n', N)
fprintf('Number of Passes           %d passes \n', Np)
fprintf('Number of Tubes/Pass       %.3f tubes/pass \n', Npass)
fprintf('Condenser Duty             %.3f kW \n', Q)
fprintf('Mass Flow Rate of CW       %.3f kg/s \n', mcw)
fprintf('Tfoi Guess                 %.3f \n', Tfoi)
fprintf('Tfoo Guess                 %.3f \n', Tfoo)