#ifndef WINDOW_
#define WINDOW_
// window.c
// this file was generated automatically using function create_signals.m
// Created on 09-Jul-2014 13:54:55

#define LEN 1024
#pragma DATA_SECTION(hanning,".EXT_RAM")
float hanning[1024] = {0.00000000e+00,1.83835134e-08,7.35333572e-08,1.65447460e-07,2.94122373e-07,4.59553235e-07,
6.61733793e-07,9.00656516e-07,1.17631225e-06,1.48869083e-06,1.83778036e-06,2.22356766e-06,
2.64603841e-06,3.10517657e-06,3.60096465e-06,4.13338421e-06,4.70241548e-06,5.30803618e-06,
5.95022448e-06,6.62895582e-06,7.34420428e-06,8.09594349e-06,8.88414525e-06,9.70877863e-06,
1.05698145e-05,1.14672193e-05,1.24009584e-05,1.33709991e-05,1.43773023e-05,1.54198315e-05,
1.64985468e-05,1.76134090e-05,1.87643745e-05,1.99514016e-05,2.11744427e-05,2.24334544e-05,
2.37283894e-05,2.50591966e-05,2.64258288e-05,2.78282314e-05,2.92663535e-05,3.07401424e-05,
3.22495398e-05,3.37944912e-05,3.53749347e-05,3.69908157e-05,3.86420725e-05,4.03286394e-05,
4.20504548e-05,4.38074530e-05,4.55995687e-05,4.74267363e-05,4.92888867e-05,5.11859471e-05,
5.31178448e-05,5.50845143e-05,5.70858756e-05,5.91218522e-05,6.11923679e-05,6.32973533e-05,
6.54367177e-05,6.76103882e-05,6.98182776e-05,7.20603130e-05,7.43363926e-05,7.66464436e-05,
7.89903788e-05,8.13681036e-05,8.37795378e-05,8.62245797e-05,8.87031420e-05,9.12151299e-05,
9.37604564e-05,9.63390194e-05,9.89507171e-05,1.01595462e-04,1.04273146e-04,1.06983665e-04,
1.09726934e-04,1.12502836e-04,1.15311268e-04,1.18152122e-04,1.21025296e-04,1.23930688e-04,
1.26868166e-04,1.29837645e-04,1.32839006e-04,1.35872135e-04,1.38936899e-04,1.42033212e-04,
1.45160942e-04,1.48319974e-04,1.51510190e-04,1.54731460e-04,1.57983683e-04,1.61266726e-04,
1.64580459e-04,1.67924765e-04,1.71299514e-04,1.74704575e-04,1.78139831e-04,1.81605152e-04,
1.85100391e-04,1.88625432e-04,1.92180145e-04,1.95764384e-04,1.99378017e-04,2.03020914e-04,
2.06692945e-04,2.10393948e-04,2.14123807e-04,2.17882363e-04,2.21669485e-04,2.25485026e-04,
2.29328842e-04,2.33200801e-04,2.37100743e-04,2.41028523e-04,2.44983996e-04,2.48967001e-04,
2.52977421e-04,2.57015083e-04,2.61079811e-04,2.65171489e-04,2.69289972e-04,2.73435056e-04,
2.77606625e-04,2.81804503e-04,2.86028546e-04,2.90278578e-04,2.94554455e-04,2.98856001e-04,
3.03183071e-04,3.07535491e-04,3.11913085e-04,3.16315709e-04,3.20743216e-04,3.25195375e-04,
3.29672097e-04,3.34173150e-04,3.38698388e-04,3.43247666e-04,3.47820751e-04,3.52417526e-04,
3.57037818e-04,3.61681421e-04,3.66348162e-04,3.71037866e-04,3.75750387e-04,3.80485551e-04,
3.85243125e-04,3.90022964e-04,3.94824892e-04,3.99648736e-04,4.04494291e-04,4.09361382e-04,
4.14249807e-04,4.19159420e-04,4.24090016e-04,4.29041422e-04,4.34013433e-04,4.39005875e-04,
4.44018573e-04,4.49051295e-04,4.54103894e-04,4.59176139e-04,4.64267883e-04,4.69378923e-04,
4.74509026e-04,4.79658076e-04,4.84825810e-04,4.90012055e-04,4.95216635e-04,5.00439317e-04,
5.05679927e-04,5.10938233e-04,5.16214117e-04,5.21507347e-04,5.26817690e-04,5.32144913e-04,
5.37488959e-04,5.42849477e-04,5.48226351e-04,5.53619291e-04,5.59028238e-04,5.64452843e-04,
5.69892989e-04,5.75348444e-04,5.80819033e-04,5.86304464e-04,5.91804623e-04,5.97319275e-04,
6.02848188e-04,6.08391187e-04,6.13948039e-04,6.19518571e-04,6.25102490e-04,6.30699680e-04,
6.36309909e-04,6.41932886e-04,6.47568493e-04,6.53216499e-04,6.58876670e-04,6.64548774e-04,
6.70232694e-04,6.75928080e-04,6.81634818e-04,6.87352673e-04,6.93081354e-04,6.98820746e-04,
7.04570615e-04,7.10330729e-04,7.16100796e-04,7.21880759e-04,7.27670267e-04,7.33469147e-04,
7.39277166e-04,7.45094148e-04,7.50919804e-04,7.56754016e-04,7.62596435e-04,7.68446946e-04,
7.74305314e-04,7.80171249e-04,7.86044635e-04,7.91925122e-04,7.97812594e-04,8.03706818e-04,
8.09607503e-04,8.15514533e-04,8.21427617e-04,8.27346521e-04,8.33271013e-04,8.39200919e-04,
8.45136005e-04,8.51076038e-04,8.57020845e-04,8.62970075e-04,8.68923613e-04,8.74881225e-04,
8.80842621e-04,8.86807684e-04,8.92776065e-04,8.98747647e-04,9.04722081e-04,9.10699309e-04,
9.16678982e-04,9.22660867e-04,9.28644848e-04,9.34630632e-04,9.40617931e-04,9.46606626e-04,
9.52596485e-04,9.58587159e-04,9.64578590e-04,9.70570429e-04,9.76562500e-04,9.82554629e-04,
9.88546410e-04,9.94537841e-04,1.00052857e-03,1.00651837e-03,1.01250701e-03,1.01849437e-03,
1.02448009e-03,1.03046407e-03,1.03644608e-03,1.04242575e-03,1.04840286e-03,1.05437741e-03,
1.06034894e-03,1.06631732e-03,1.07228232e-03,1.07824383e-03,1.08420139e-03,1.09015487e-03,
1.09610416e-03,1.10204890e-03,1.10798900e-03,1.11392408e-03,1.11985393e-03,1.12577854e-03,
1.13169744e-03,1.13761052e-03,1.14351744e-03,1.14941818e-03,1.15531241e-03,1.16119988e-03,
1.16708037e-03,1.17295375e-03,1.17881969e-03,1.18467805e-03,1.19052851e-03,1.19637104e-03,
1.20220520e-03,1.20803085e-03,1.21384778e-03,1.21965585e-03,1.22545473e-03,1.23124430e-03,
1.23702420e-03,1.24279433e-03,1.24855433e-03,1.25430420e-03,1.26004359e-03,1.26577239e-03,
1.27149024e-03,1.27719692e-03,1.28289231e-03,1.28857617e-03,1.29424839e-03,1.29990850e-03,
1.30555651e-03,1.31119206e-03,1.31681515e-03,1.32242532e-03,1.32802245e-03,1.33360643e-03,
1.33917690e-03,1.34473376e-03,1.35027675e-03,1.35580567e-03,1.36132038e-03,1.36682054e-03,
1.37230603e-03,1.37777650e-03,1.38323195e-03,1.38867216e-03,1.39409676e-03,1.39950565e-03,
1.40489871e-03,1.41027558e-03,1.41563604e-03,1.42098009e-03,1.42630737e-03,1.43161765e-03,
1.43691082e-03,1.44218677e-03,1.44744513e-03,1.45268568e-03,1.45790842e-03,1.46311300e-03,
1.46829919e-03,1.47346698e-03,1.47861592e-03,1.48374611e-03,1.48885709e-03,1.49394886e-03,
1.49902108e-03,1.50407373e-03,1.50910649e-03,1.51411910e-03,1.51911157e-03,1.52408355e-03,
1.52903493e-03,1.53396558e-03,1.53887516e-03,1.54376368e-03,1.54863077e-03,1.55347632e-03,
1.55830011e-03,1.56310201e-03,1.56788190e-03,1.57263945e-03,1.57737464e-03,1.58208713e-03,
1.58677681e-03,1.59144355e-03,1.59608724e-03,1.60070742e-03,1.60530419e-03,1.60987733e-03,
1.61442661e-03,1.61895179e-03,1.62345287e-03,1.62792963e-03,1.63238181e-03,1.63680932e-03,
1.64121191e-03,1.64558948e-03,1.64994190e-03,1.65426906e-03,1.65857060e-03,1.66284642e-03,
1.66709651e-03,1.67132053e-03,1.67551835e-03,1.67968997e-03,1.68383506e-03,1.68795348e-03,
1.69204513e-03,1.69610989e-03,1.70014764e-03,1.70415803e-03,1.70814106e-03,1.71209651e-03,
1.71602424e-03,1.71992416e-03,1.72379613e-03,1.72763993e-03,1.73145556e-03,1.73524267e-03,
1.73900125e-03,1.74273108e-03,1.74643204e-03,1.75010413e-03,1.75374700e-03,1.75736065e-03,
1.76094484e-03,1.76449958e-03,1.76802464e-03,1.77151989e-03,1.77498523e-03,1.77842041e-03,
1.78182544e-03,1.78520021e-03,1.78854458e-03,1.79185823e-03,1.79514126e-03,1.79839355e-03,
1.80161477e-03,1.80480501e-03,1.80796406e-03,1.81109179e-03,1.81418809e-03,1.81725284e-03,
1.82028604e-03,1.82328734e-03,1.82625686e-03,1.82919437e-03,1.83209975e-03,1.83497288e-03,
1.83781376e-03,1.84062216e-03,1.84339809e-03,1.84614130e-03,1.84885191e-03,1.85152958e-03,
1.85417430e-03,1.85678597e-03,1.85936457e-03,1.86190987e-03,1.86442188e-03,1.86690048e-03,
1.86934543e-03,1.87175686e-03,1.87413464e-03,1.87647855e-03,1.87878858e-03,1.88106473e-03,
1.88330677e-03,1.88551459e-03,1.88768830e-03,1.88982766e-03,1.89193268e-03,1.89400313e-03,
1.89603912e-03,1.89804053e-03,1.90000713e-03,1.90193905e-03,1.90383615e-03,1.90569821e-03,
1.90752547e-03,1.90931757e-03,1.91107451e-03,1.91279640e-03,1.91448291e-03,1.91613415e-03,
1.91775011e-03,1.91933056e-03,1.92087551e-03,1.92238484e-03,1.92385865e-03,1.92529673e-03,
1.92669919e-03,1.92806579e-03,1.92939665e-03,1.93069153e-03,1.93195057e-03,1.93317363e-03,
1.93436060e-03,1.93551159e-03,1.93662650e-03,1.93770521e-03,1.93874771e-03,1.93975400e-03,
1.94072409e-03,1.94165774e-03,1.94255519e-03,1.94341619e-03,1.94424088e-03,1.94502901e-03,
1.94578082e-03,1.94649608e-03,1.94717478e-03,1.94781693e-03,1.94842264e-03,1.94899156e-03,
1.94952404e-03,1.95001985e-03,1.95047900e-03,1.95090147e-03,1.95128727e-03,1.95163628e-03,
1.95194874e-03,1.95222429e-03,1.95246330e-03,1.95266539e-03,1.95283093e-03,1.95295957e-03,
1.95305143e-03,1.95310661e-03,1.95312500e-03,1.95310661e-03,1.95305143e-03,1.95295957e-03,
1.95283093e-03,1.95266539e-03,1.95246330e-03,1.95222429e-03,1.95194874e-03,1.95163628e-03,
1.95128727e-03,1.95090147e-03,1.95047900e-03,1.95001985e-03,1.94952404e-03,1.94899156e-03,
1.94842264e-03,1.94781693e-03,1.94717478e-03,1.94649608e-03,1.94578082e-03,1.94502901e-03,
1.94424088e-03,1.94341619e-03,1.94255519e-03,1.94165774e-03,1.94072409e-03,1.93975400e-03,
1.93874771e-03,1.93770521e-03,1.93662650e-03,1.93551159e-03,1.93436060e-03,1.93317363e-03,
1.93195057e-03,1.93069153e-03,1.92939665e-03,1.92806579e-03,1.92669919e-03,1.92529673e-03,
1.92385865e-03,1.92238484e-03,1.92087551e-03,1.91933056e-03,1.91775011e-03,1.91613415e-03,
1.91448291e-03,1.91279640e-03,1.91107451e-03,1.90931757e-03,1.90752547e-03,1.90569821e-03,
1.90383615e-03,1.90193905e-03,1.90000713e-03,1.89804053e-03,1.89603912e-03,1.89400313e-03,
1.89193268e-03,1.88982766e-03,1.88768830e-03,1.88551459e-03,1.88330677e-03,1.88106473e-03,
1.87878858e-03,1.87647855e-03,1.87413464e-03,1.87175686e-03,1.86934543e-03,1.86690048e-03,
1.86442188e-03,1.86190987e-03,1.85936457e-03,1.85678597e-03,1.85417430e-03,1.85152958e-03,
1.84885191e-03,1.84614130e-03,1.84339809e-03,1.84062216e-03,1.83781376e-03,1.83497288e-03,
1.83209975e-03,1.82919437e-03,1.82625686e-03,1.82328734e-03,1.82028604e-03,1.81725284e-03,
1.81418809e-03,1.81109179e-03,1.80796406e-03,1.80480501e-03,1.80161477e-03,1.79839355e-03,
1.79514126e-03,1.79185823e-03,1.78854458e-03,1.78520021e-03,1.78182544e-03,1.77842041e-03,
1.77498523e-03,1.77151989e-03,1.76802464e-03,1.76449958e-03,1.76094484e-03,1.75736065e-03,
1.75374700e-03,1.75010413e-03,1.74643204e-03,1.74273108e-03,1.73900125e-03,1.73524267e-03,
1.73145556e-03,1.72763993e-03,1.72379613e-03,1.71992416e-03,1.71602424e-03,1.71209651e-03,
1.70814106e-03,1.70415803e-03,1.70014764e-03,1.69610989e-03,1.69204513e-03,1.68795348e-03,
1.68383506e-03,1.67968997e-03,1.67551835e-03,1.67132053e-03,1.66709651e-03,1.66284642e-03,
1.65857060e-03,1.65426906e-03,1.64994190e-03,1.64558948e-03,1.64121191e-03,1.63680932e-03,
1.63238181e-03,1.62792963e-03,1.62345287e-03,1.61895179e-03,1.61442661e-03,1.60987733e-03,
1.60530419e-03,1.60070742e-03,1.59608724e-03,1.59144355e-03,1.58677681e-03,1.58208713e-03,
1.57737464e-03,1.57263945e-03,1.56788190e-03,1.56310201e-03,1.55830011e-03,1.55347632e-03,
1.54863077e-03,1.54376368e-03,1.53887516e-03,1.53396558e-03,1.52903493e-03,1.52408355e-03,
1.51911157e-03,1.51411910e-03,1.50910649e-03,1.50407373e-03,1.49902108e-03,1.49394886e-03,
1.48885709e-03,1.48374611e-03,1.47861592e-03,1.47346698e-03,1.46829919e-03,1.46311300e-03,
1.45790842e-03,1.45268568e-03,1.44744513e-03,1.44218677e-03,1.43691082e-03,1.43161765e-03,
1.42630737e-03,1.42098009e-03,1.41563604e-03,1.41027558e-03,1.40489871e-03,1.39950565e-03,
1.39409676e-03,1.38867216e-03,1.38323195e-03,1.37777650e-03,1.37230603e-03,1.36682054e-03,
1.36132038e-03,1.35580567e-03,1.35027675e-03,1.34473376e-03,1.33917690e-03,1.33360643e-03,
1.32802245e-03,1.32242532e-03,1.31681515e-03,1.31119206e-03,1.30555651e-03,1.29990850e-03,
1.29424839e-03,1.28857617e-03,1.28289231e-03,1.27719692e-03,1.27149024e-03,1.26577239e-03,
1.26004359e-03,1.25430420e-03,1.24855433e-03,1.24279433e-03,1.23702420e-03,1.23124430e-03,
1.22545473e-03,1.21965585e-03,1.21384778e-03,1.20803085e-03,1.20220520e-03,1.19637104e-03,
1.19052851e-03,1.18467805e-03,1.17881969e-03,1.17295375e-03,1.16708037e-03,1.16119988e-03,
1.15531241e-03,1.14941818e-03,1.14351744e-03,1.13761052e-03,1.13169744e-03,1.12577854e-03,
1.11985393e-03,1.11392408e-03,1.10798900e-03,1.10204890e-03,1.09610416e-03,1.09015487e-03,
1.08420139e-03,1.07824383e-03,1.07228232e-03,1.06631732e-03,1.06034894e-03,1.05437741e-03,
1.04840286e-03,1.04242575e-03,1.03644608e-03,1.03046407e-03,1.02448009e-03,1.01849437e-03,
1.01250701e-03,1.00651837e-03,1.00052857e-03,9.94537841e-04,9.88546410e-04,9.82554629e-04,
9.76562500e-04,9.70570429e-04,9.64578590e-04,9.58587159e-04,9.52596485e-04,9.46606626e-04,
9.40617931e-04,9.34630632e-04,9.28644848e-04,9.22660867e-04,9.16678982e-04,9.10699309e-04,
9.04722081e-04,8.98747647e-04,8.92776065e-04,8.86807684e-04,8.80842621e-04,8.74881225e-04,
8.68923613e-04,8.62970075e-04,8.57020845e-04,8.51076038e-04,8.45136005e-04,8.39200919e-04,
8.33271013e-04,8.27346521e-04,8.21427617e-04,8.15514533e-04,8.09607503e-04,8.03706818e-04,
7.97812594e-04,7.91925122e-04,7.86044635e-04,7.80171249e-04,7.74305314e-04,7.68446946e-04,
7.62596435e-04,7.56754016e-04,7.50919804e-04,7.45094148e-04,7.39277166e-04,7.33469147e-04,
7.27670267e-04,7.21880759e-04,7.16100796e-04,7.10330729e-04,7.04570615e-04,6.98820746e-04,
6.93081354e-04,6.87352673e-04,6.81634818e-04,6.75928080e-04,6.70232694e-04,6.64548774e-04,
6.58876670e-04,6.53216499e-04,6.47568493e-04,6.41932886e-04,6.36309909e-04,6.30699680e-04,
6.25102490e-04,6.19518571e-04,6.13948039e-04,6.08391187e-04,6.02848188e-04,5.97319275e-04,
5.91804623e-04,5.86304464e-04,5.80819033e-04,5.75348444e-04,5.69892989e-04,5.64452843e-04,
5.59028238e-04,5.53619291e-04,5.48226351e-04,5.42849477e-04,5.37488959e-04,5.32144913e-04,
5.26817690e-04,5.21507347e-04,5.16214117e-04,5.10938233e-04,5.05679927e-04,5.00439317e-04,
4.95216635e-04,4.90012055e-04,4.84825810e-04,4.79658076e-04,4.74509026e-04,4.69378923e-04,
4.64267883e-04,4.59176139e-04,4.54103894e-04,4.49051295e-04,4.44018573e-04,4.39005875e-04,
4.34013433e-04,4.29041422e-04,4.24090016e-04,4.19159420e-04,4.14249807e-04,4.09361382e-04,
4.04494291e-04,3.99648736e-04,3.94824892e-04,3.90022964e-04,3.85243125e-04,3.80485551e-04,
3.75750387e-04,3.71037866e-04,3.66348162e-04,3.61681421e-04,3.57037818e-04,3.52417526e-04,
3.47820751e-04,3.43247666e-04,3.38698388e-04,3.34173150e-04,3.29672097e-04,3.25195375e-04,
3.20743216e-04,3.16315709e-04,3.11913085e-04,3.07535491e-04,3.03183071e-04,2.98856001e-04,
2.94554455e-04,2.90278578e-04,2.86028546e-04,2.81804503e-04,2.77606625e-04,2.73435056e-04,
2.69289972e-04,2.65171489e-04,2.61079811e-04,2.57015083e-04,2.52977421e-04,2.48967001e-04,
2.44983996e-04,2.41028523e-04,2.37100743e-04,2.33200801e-04,2.29328842e-04,2.25485026e-04,
2.21669485e-04,2.17882363e-04,2.14123807e-04,2.10393948e-04,2.06692945e-04,2.03020914e-04,
1.99378017e-04,1.95764384e-04,1.92180145e-04,1.88625432e-04,1.85100391e-04,1.81605152e-04,
1.78139831e-04,1.74704575e-04,1.71299514e-04,1.67924765e-04,1.64580459e-04,1.61266726e-04,
1.57983683e-04,1.54731460e-04,1.51510190e-04,1.48319974e-04,1.45160942e-04,1.42033212e-04,
1.38936899e-04,1.35872135e-04,1.32839006e-04,1.29837645e-04,1.26868166e-04,1.23930688e-04,
1.21025296e-04,1.18152122e-04,1.15311268e-04,1.12502836e-04,1.09726934e-04,1.06983665e-04,
1.04273146e-04,1.01595462e-04,9.89507171e-05,9.63390194e-05,9.37604564e-05,9.12151299e-05,
8.87031420e-05,8.62245797e-05,8.37795378e-05,8.13681036e-05,7.89903788e-05,7.66464436e-05,
7.43363926e-05,7.20603130e-05,6.98182776e-05,6.76103882e-05,6.54367177e-05,6.32973533e-05,
6.11923679e-05,5.91218522e-05,5.70858756e-05,5.50845143e-05,5.31178448e-05,5.11859471e-05,
4.92888867e-05,4.74267363e-05,4.55995687e-05,4.38074530e-05,4.20504548e-05,4.03286394e-05,
3.86420725e-05,3.69908157e-05,3.53749347e-05,3.37944912e-05,3.22495398e-05,3.07401424e-05,
2.92663535e-05,2.78282314e-05,2.64258288e-05,2.50591966e-05,2.37283894e-05,2.24334544e-05,
2.11744427e-05,1.99514016e-05,1.87643745e-05,1.76134090e-05,1.64985468e-05,1.54198315e-05,
1.43773023e-05,1.33709991e-05,1.24009584e-05,1.14672193e-05,1.05698145e-05,9.70877863e-06,
8.88414525e-06,8.09594349e-06,7.34420428e-06,6.62895582e-06,5.95022448e-06,5.30803618e-06,
4.70241548e-06,4.13338421e-06,3.60096465e-06,3.10517657e-06,2.64603841e-06,2.22356766e-06,
1.83778036e-06,1.48869083e-06,1.17631225e-06,9.00656516e-07,6.61733793e-07,4.59553235e-07,
2.94122373e-07,1.65447460e-07,7.35333572e-08,1.83835134e-08

};


#endif /*WINDOW_*/