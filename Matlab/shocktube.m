 function [x_mesh ,u,a,rho ,T,p] = shocktube(time ,p1 ,p4 ,T1 ,T4)
 % SHOCKTUBE Analytical solution for unsteady wave motion in a shock tube .
 % [ X_MESH ,U,A,RHO ,T,P] = SHOCKTUBE (TIME ,P1 ,P4 ,T1 ,T4) solves the shock
 % tube problem analytically . The diaphragm is placed at 15.24 cm ,
 % pressure and temperature to the left of it are P4 and T4 , to the right
 % of it P1 and T1 (at time zero ). The function returns X_MESH (1000
 % equally spaced points between 0 and 30.48 cm) and the mass velocity U,
 % the local speed of sound A, the density RHO , the temperature T and the
 % pressure P at time TIME for further analysis .
 %
 % If only one argument is given ( TIME ), P1 , P4 , T1 and T4 are set to
 % default values .
 %
 % SHOCKTUBE only treats right running shock waves and left running
 % expansion waves .
 %
 % The gas in the tube is air with corresponding R and gamma ; it is
 % treated as inviscid .

 % *************************************************************************
 % This m- file is part of the Master Thesis
 %
 % " Simulation and validation of compressible flow in nozzle geometries and
 % validation of OpenFOAM for this application "
 %
 % by Benjamin Wuethrich , MSc student of Computational Science and
 % Engineering at ETH Zurich .
 %
 % Work carried out at ABB Corporate Research in Baden - Daettwil from
 % 15/04/07 until 14/09/07.
 %
 % Contact : benjamin . wuethrich@alumni . ethz .ch
 % *************************************************************************

 % Parse arguments .
 if nargin ==1
 % Set defaults .
 p1 = 6.897e3; % Lower pressure ( right chamber ) in [Pa ].
 p4 = 6.897e4; % Higher pressure ( left chamber ) in [Pa ].
 T1 = 231.11; % Lower temperature ( right chamber ) in [K].
 T4 = 288.89; % Higher temperature ( left chamber ) in [K].
 elseif ( nargin ==5) && (p1 > p4)
 error ('This would be a left - running shock wave , which is not supported.')
 elseif nargin ==5
 % Everything fine .
 else
 error ('Wrong number of arguments specified.')
 end

 % Set constants .
 gamma = 1.4; % Heat capacity ratio of air.
 R = 287.05; % Specific gas constant of air in [J/( kg*K )].
 L1 = 0.1524; % Initial position of the diaphragm in [m].$

 % Calculate speeds of sound and densities .
 a1 = sqrt ( gamma *R*T1 );
 a4 = sqrt ( gamma *R*T4 );
 rho1 = p1 /(R*T1 );
 rho4 = p4 /(R*T4 );

 % Calculate p2/p1 , get p2.
 p2p1 = fzero (@( p2p1 ) p2p1 * (1 - (( gamma -1)*( a1/a4 )*( p2p1 -1)) / ...
 sqrt (2* gamma *(2* gamma + ( gamma +1)*( p2p1 -1))) )^( -2* gamma /( gamma -1))...
 - p4/p1 ,( p4/p1 )/2);
 p2 = p2p1 * p1;

 % Calculate incident shock properties .
 T2 = T1 * p2/p1 * ((( gamma +1)/( gamma -1) + p2/p1) /...
 (1 + ( gamma +1)/( gamma -1)* p2/p1 ));
 rho2 = rho1 * (1 + ( gamma +1)/( gamma -1)* p2/p1) /...
 (( gamma +1)/( gamma -1) + p2/p1 );
 a2 = sqrt ( gamma *R*T2 );
 % Wave velocity of moving shock .
 W = a1 * sqrt (( gamma +1)/(2* gamma ) * (p2/p1 - 1) + 1);
 % Mass motion behind the wave .
 u_p = a1/ gamma * (p2/p1 - 1) * sqrt ((2* gamma /( gamma +1)) /...
 (p2/p1 + (gamma -1)/( gamma +1)));

 % Pressure and velocity to the right of the expansion wave ( constant
 % across the contact surface ).
 p3 = p2;
 u2 = u_p;
 u3 = u_p;

 % Other properties behind the expansion wave .
 rho3 = rho4 * (p3/p4 )^(1/ gamma );
 T3 = T4 * (p3/p4 )^(( gamma -1)/ gamma );
 a3 = sqrt ( gamma *R*T3 );

 % Define mesh .
 x_mesh = linspace (0 ,2* L1 ,100);

 % Initialise vectors for all the quantities .
 u = zeros ( size ( x_mesh ));
 a = zeros ( size ( x_mesh ));
 rho = zeros ( size ( x_mesh ));
 T = zeros ( size ( x_mesh ));
 p = zeros ( size ( x_mesh ));

 % Calculate boundaries of different zones .
 % Boundary between leftmost driver gas and expansion wave .
 x4_exp = L1 - time *a4;
 % Boundary between expansion wave and lower pressure driver gas.
 exp_x3 = L1 + time *( u3 - a3 );
 % Boundary between driver gas and driven gas.
 x3_x2 = L1 + time *u_p;
 % Location of the shock wave .
 x2_x1 = L1 + time *W;

 % Iterate through x_mesh and fill in all the quantities .
 for i = 1: length ( x_mesh )
 if x_mesh (i) < x4_exp
 % We are in region 4.
 u(i) = 0;
 a(i) = a4;
 rho(i) = rho4 ;
 T(i) = T4;
 p(i) = p4;
 elseif x_mesh (i) < exp_x3
 % We are in the expansion wave .
 [u(i),a(i),rho(i),T(i),p(i)] = expansion_wave ( x_mesh (i)-L1 );
 elseif x_mesh (i) < x3_x2
 % We are in region 3.
 u(i) = u3;
 a(i) = a3;
 rho(i) = rho3 ;
 T(i) = T3;
 p(i) = p3;
 elseif x_mesh (i) < x2_x1
 % We are in region 2.
 u(i) = u2;
 a(i) = a2;
 rho(i) = rho2 ;
 T(i) = T2;
 p(i) = p2;
 else
 % We are in region 1.
 u(i) = 0;
 a(i) = a1;
 rho(i) = rho1 ;
 T(i) = T1;
 p(i) = p1;
 end
 end

 function [u_exp ,a_exp , rho_exp ,T_exp , p_exp ] = expansion_wave (x)
 % Calculate properties within expansion wave . The expressions here
 % are valid for -a4 <= x/ time <= u3 - a3.
 % Calculate mass velocity .
 u_exp = 2/( gamma +1) * (a4 + x/ time );
 % Calculate speed of sound .
 a_exp = a4 * (1 - (gamma -1)/2 * u_exp /a4 );
 % Calculate temperature .
 T_exp = T4 * (1 - (gamma -1)/2 * u_exp /a4 )^2;
 % Calculate pressure .
 p_exp = p4 * (1 - (gamma -1)/2 * u_exp /a4 )^(2* gamma /( gamma -1));
 % Calculate density .
 rho_exp = rho4 * (1 - (gamma -1)/2 * u_exp /a4 )^(2/( gamma -1));
 end
 end