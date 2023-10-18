clc
clear
close all

%Case I

r1 = 9;
b1 = 1;
r2 = 9;
b2 = 1;

[prior1, likeR1, likeB1, postR1, postB1] =  calculatePosteriors(r1, b1, r2, b2);

mlDecVec1 = genMlDecVec(likeR1, likeB1);
mapDecVec1 = genMapDecVec(postR1, postB1);

Perr_ml1 = genMlErr(mlDecVec1, likeB1, likeR1, prior1);
Perr_map1 = genMapErr(mapDecVec1, likeB1, likeR1, prior1);

[percentMlErr1, percentMapErr1] = runExpmt(r1, b1, r2, b2, mlDecVec1, mapDecVec1);



%Case II

r1 = 6;
b1 = 4;
r2 = 9;
b2 = 1;

[prior2, likeR2, likeB2, postR2, postB2] =  calculatePosteriors(r1, b1, r2, b2);

mlDecVec2 = genMlDecVec(likeR2, likeB2);
mapDecVec2 = genMapDecVec(postR2, postB2);

Perr_ml2 = genMlErr(mlDecVec2, likeB2, likeR2, prior2);
Perr_map2 = genMapErr(mapDecVec2, likeB2, likeR2, prior2);

[percentMlErr2, percentMapErr2] = runExpmt(r1, b1, r2, b2, mlDecVec2, mapDecVec2);

%Case III

r1 = 9;
b1 = 1;
r2 = 6;
b2 = 4;

[prior3, likeR3, likeB3, postR3, postB3] =  calculatePosteriors(r1, b1, r2, b2);

mlDecVec3 = genMlDecVec(likeR3, likeB3);
mapDecVec3 = genMapDecVec(postR3, postB3);

Perr_ml3 = genMlErr(mlDecVec3, likeB3, likeR3, prior3);
Perr_map3 = genMapErr(mapDecVec3, likeB3, likeR3, prior3);

[percentMlErr3, percentMapErr3] = runExpmt(r1, b1, r2, b2, mlDecVec3, mapDecVec3);


%part E

% The ML decision should NOT depend on the distribution 
% of balls in Urn I. This IS consistent because the ML
% decision only depends on the likelihoods.
% The likelihood calculations do not depend on the the balls
% in URN I, so neither does the ML decision

%part E

% The MAP decision SHOULD depend on the distribution 
% of balls in Urn I because it depends on the posteriors
% The posterior calculations factor in the prior functions
% and the priors depend on the balls in URN I, 
% so therefore the MAP decision also depends on the URN I balls

%part G

% The ML errors are almost exact!
% See percentMlErr(1-3) [experimental] and Perr_ml(1-3) [theory]

% The MAP errors are wrong!
% I definetly calculated it wrong
% See percentMapErr(1-3) [experimental] and Perr_map(1-3) [theory]

function [percentMlErr, percentMapErr] = runExpmt(r1, b1, r2, b2, mlDecVec, mapDecVec)

    mlErrs = 0;

    mapErrs = 0;

    for N = 1:10000
    
        urnI = [ones(1, b1), 2 * ones(1, r1)];
    
        cUrnI = urnI(randi(length(urnI)));
    
        urnII = [ones(1, b2 + (cUrnI == 1)), 2 * ones(1, r2 + (cUrnI == 2))];
    
        cUrnII = urnII(randi(length(urnII)));
    
        if cUrnII == 1
    
            if mlDecVec(1) ~= cUrnI
            
                mlErrs = mlErrs + 1;
            
            end
        
            if mapDecVec(1) ~= cUrnI
            
                mapErrs = mapErrs + 1;
            
            end
        
        end
    
        if cUrnII == 2
    
            if mlDecVec(2) ~= cUrnI
            
                mlErrs = mlErrs + 1;
            
            end
        
            if mapDecVec(2) ~= cUrnI
            
                mapErrs = mapErrs + 1;
            
            end  
        end  
    end

    percentMlErr = mlErrs/N;
    percentMapErr = mapErrs/N;
    
end

function P_err_ml = genMlErr(mlDecVec, likeB, likeR, prior)

    P_err_ml1 = 0;
    P_err_ml2 = 0;
    
    %we get blue from urn II

    if mlDecVec(1) == 1 %chose blue
        
        P_err_ml1 = (1 - likeB(1))*prior(1);
        
    end
    
    if mlDecVec(1) == 2 %chose red
        
        P_err_ml1 = (1 - likeB(2))*prior(2);
        
    end
    
    %we get rd from urn II
    
    if mlDecVec(2) == 1 %chose blue
    
        P_err_ml2 = (1 - likeR(1))*prior(1);
    
    end

    if mlDecVec(2) == 2 %chose red
    
        P_err_ml2 = (1 - likeR(2))*prior(2);
    
    end
    
    P_err_ml = P_err_ml1 + P_err_ml2;

end

function P_err_map = genMapErr(mapDecVec, likeB, likeR, prior)

    P_err_map1 = 0;
    P_err_map2 = 0;
    
    %we get blue from urn II

    if mapDecVec(1) == 1 %chose blue
        
        P_err_map1 = (1-likeB(1))*prior(1);
        
    end
    
    if mapDecVec(1) == 2 %chose red
        
        P_err_map1 = (1-likeB(2))*prior(2);
        
    end
    
    %we get rd from urn II
    
    if mapDecVec(2) == 1 %chose blue
    
        P_err_map2 = (1-likeR(1))*prior(1);
    
    end

    if mapDecVec(2) == 2 %chose red
    
        P_err_map2 = (1-likeR(2))*prior(2);
    
    end
    
    P_err_map = 1 - (P_err_map1 + P_err_map2);

end

function mldv = genMlDecVec(likeR, likeB)

    mldv = zeros(1, 2);
    
    [maxR, iR] = max(likeR);
    [maxB, iB] = max(likeB);
    
    mldv(1) = iB;
    mldv(2) = iR;  

end


function mapdv = genMapDecVec(postR, postB)

    mapdv = zeros(1, 2);
    
    [maxR, iR] = max(postR);
    [maxB, iB] = max(postB);
    
    mapdv(1) = iB;
    mapdv(2) = iR;    

end

function [prior, likeR, likeB, postR, postB] = calculatePosteriors(r1, b1, r2, b2)

    %prior models

    pRi = r1/(r1+b1);
    pBi = b1/(r1+b1);

    %likelihood functions

    pRii_gRi = (r2+1)/(r2+b2+1);
    pRii_gBi = (r2)/(r2+b2+1);
    pBii_gRi = (b2)/(r2+b2+1);
    pBii_gBi = (b2+1)/(r2+b2+1);

    %posterior distributions

    pRi_gRii = (pRii_gRi * pRi) / ((pRii_gRi * pRi) + (pRii_gBi * pBi));
    pRi_gBii = (pBii_gRi * pRi) / ((pBii_gRi * pRi) + (pBii_gBi * pBi));

    pBi_gRii = (pRii_gBi * pBi) / ((pRii_gRi * pRi) + (pRii_gBi * pBi));
    pBi_gBii = (pBii_gBi * pBi) / ((pBii_gRi * pRi) + (pBii_gBi * pBi));
    
    prior = [pBi, pRi];
    
    likeR = [pRii_gBi, pRii_gRi];
    likeB = [pBii_gBi, pBii_gRi];
    
    postR = [pBi_gRii, pRi_gRii];
    postB = [pBi_gBii, pRi_gBii];
    
end