%//CFWH #3 Sizing//

%input values below:
Nexch = 1;
ms = 55.56*0.068; %steam mass flow rate of CFWH
Rfi = 0.00015; Rfo = 0.00013; %fouling factors
Tho = 151.836; Thi = 505.39; Tmo = 148.836 ; Tmi = 102.97; %temperature values
hin = 3496.02; hout = 640.185; %enthalpy values entering and exiting the CFWH at states 2 and 21
u_cw = 2; %velocity of cooling water (m = pAv)

%Cooling water properties T = 125.9C @ 10MPa

k = 0.687; v = 2.35e-7; Prw = 1.436; p = 943.24; Cpw = 4.23; u = 2.33e-4;

%Steam properties T = 263.943C @ 5MPa

hf = 639.98; hfg = 2108.64; Prs = 1.145; pl = 587.88; Cps = 3.096; ul = 1.76e-4; kl = 0.659; 

%properties of tubes (1-1/4in Schedule 40 square pitch):
 
Di = 0.035052; Di_inches = 1.38; Do = 0.04216; Do_inches = 1.660; tw = 0.003556; tw_inches = 0.140; PT = 0.03969;  
Np = 6; %Number of passes in the shell
kw = 48.9; %thermal conductivity of 1 Cr–V(0.2% C, 1.02% Cr, 0.15% V) alloy steel

%Start code:

Tlm = ((Tho - Tmo) - (Thi - Tmi)) / log((Tho - Tmo) / (Thi - Tmi));

Tm = (Tmi + Tmo)/2;

Q = ms*(hin - hout); %required rate of heat transfer from steam to cooling water in kW

mcw = Q / (Cpw*(Tmo -Tmi));
mcw = 55.56; %mass flow rate of the cooling water in kg/s

Ntube = 4*mcw / (u_cw*p*pi*Di^2); %number of tubes/pass
Ntube = round(Ntube, 0);
Nr = Ntube*Np; %number of tubes required

N = 1008; Npass = 1008/6; Ds = 1.524; Ds_inches = 60; %shell tube layouts TEMA standard

ucw = 4*mcw / (p*pi*Di^2*Npass); %velocity of cooling water inside condenser tubes for the shell diameter
ReD = ucw*p*Di / u; %tube side Reynolds number
ReD = round(ReD, 0);

f = (0.790*log(ReD)- 1.64)^-2; %frictional factor for Nuesslet number equation below
NuD = ((f/8)*(ReD - 1000)*Prs) / (1 + ((12.7*(f/8)^0.5) * (Prs^(2/3) - 1))); %Nusselt number for a smooth pipe

hi = NuD*(k/Di); %tube side convection coefficient

NT = Ds/PT; %average number of transverse tubes
NT = round(NT, 0);

Dm = (Do - Di) / (log(Do/Di));  
Rt = Rfo + ((1+hi*Rfi)/hi)*(Do/Di) + (tw/kw)*(Do/Dm); %thermal resistance (m^2*K/W)

dTin = Thi - Tmi;
dTout = Tho - Tmo;
Tfoi = 120.082; Tfoo = 0.261; %guesses for fouled surface at the inlet and outlet of the tube (C or K)

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
Ptot = (Pt + Pr)/10^3; %total pressure drop through the tubes (kPa
Ptot_psi = Ptot*0.145038; %convert to psi

%Pumping power required for the cooling water:

PWR = ((mcw*Ptot)/(p*0.85)); %power in kW

%Display sizing values:

fprintf('Final Design Summary: \n\n\n')
fprintf('Number of Heat Exchangers  %d \n', Nexch)
fprintf('Shell Diameter             %.3f m (max 60 in) \n', Ds)
fprintf('Shell Length               %.2f m (max 35ft) \n', L)
fprintf('Inner Diameter             %.3f m \n', Di)
fprintf('Outer Diameter             %.3f m \n', Do)
fprintf('Tube Pressure Drop         %.2f kpa \n', Ptot)
fprintf('Pumping Power              %.2f kW \n', PWR)
fprintf('Number of Tubes Required   %d tubes \n', Nr)
fprintf('Transverse Tubes           %d tubes \n', NT)
fprintf('Number of Tubes            %d tubes \n', N)
fprintf('Number of Passes           %d passes \n', Np)
fprintf('Number of Tubes/Pass       %d tubes/pass \n', Npass)
fprintf('Tfoi Guess                 %.3f \n', Tfoi)
fprintf('Tfoo Guess                 %.3f \n', Tfoo)
