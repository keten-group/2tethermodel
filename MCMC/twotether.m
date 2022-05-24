% Markov Chain Monte Carlo (MCMC) Code for a two tether system 
% One has a loop and the other doesn't. Initially loop is closed.
% Both the loop and the adhesins have off-rate k = A*exp(B*f)
% A -> k0, B ~ x0. 1/A is like the zero-force lifetime. 

%INPUT
% N = Number of trials for each force value
% fvec = Applied force

%OUTPUT
% LTavg = Average lifetime for each force
% Lifetime (LTavg) vs force (fvec) curve

%EXAMPLE
%[fvec,LTavg] = twotether(10000,[.1:0.1:2])

function [fvec,LTavg] = twotether(N,fvec)
tic 
LT = zeros(N,length(fvec));
X = length(fvec);
for i = 1:N % repeat N times for each f
    for j = 1:X % iterate over force
        %j
        [LT(i,j)]=lifetime(fvec(j)); % lifetime matrix
    end
end
toc 
 
LTavg = mean(LT,1); % take average over columns, over N trials
 
figure (1) 
 
plot(fvec,LTavg,'ok')
 
function [LT]=lifetime(f)
%time step in sec
dt=1e-8;
 
% kinetic constants
 
a = 1e-3/5  ;
b = 1/2 ;


c = 100 ;
d = 0.5/2  ;

% Initially, all bonds intact, 
% Flags to check if the bond is broken or not, 0 = intact, 1 = broken
% Bl -> loop, Ba1 -> loop-side adhesin, Ba2-> other adh

Bl = 0; 
Ba1 = 0;
Ba2 = 0;
t = 0;

while 1
    t=t+dt;
    
    if Bl == 0 
        [Bl] = bondcheck(f,a,b,dt); % test the loop if not broken previously
    end
   
    if Ba1 == 1 % adh1 broken?
        [Ba2] = bondcheck(f,c,d,dt); % adh2 carries all
    elseif Ba2 == 1 %adh2 broken?
        [Ba1] = bondcheck(f,c,d,dt); %adh 1 carries all, loop irrelevant.
    else
        
        if Bl == 0  % loop intact
            [Ba1] = bondcheck(f,c,d,dt); % adh1 carries all
            [Ba2] = bondcheck(0,c,d,dt); % adh2 no load
        elseif Bl == 1 %loop broken
           % 'loop open!' % loop open flag, if needed
            [Ba1] = bondcheck(f/2,c,d,dt); % 1/2 load for adh1 
            [Ba2] = bondcheck(f/2,c,d,dt); % 1/2 load for adh2
        end
    end
    if Ba1==1
        if Ba2==1
            %f;
            LT = t;
            break
        end
    end
end
 
function [B] = bondcheck(f,A,B,dt)
 
k = A*exp(B*f);
P = k*dt ;
r=rand;
 if P>r
     B = 1; % B equals 1 if bond broken
 else
     B=0; 
 end
