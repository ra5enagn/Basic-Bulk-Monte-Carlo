function [kxn,kyn,kzn,energy] = drift(tau,iv,kx,ky,kz,fx,fy,fz,tm1,tm2,tm3,hhm,energy,alpha)
    h = 6.625e-34;
    h_cross = h/(2*pi);
    q = 1.602e-19;
    kb= 1.38064e-23;
    

    qh1 = -1*q*tau/h_cross;
    if(iv == 1)
        dkx = qh1 * fx * tm1;
        dky = qh1 * fy * tm2;
        dkz = qh1 * fz * tm3;
    elseif (iv == 2)
        dkx = qh1 * fx * tm2;
        dky = qh1 * fy * tm1;
        dkz = qh1 * fz * tm3;
    elseif (iv == 3)
        dkx = qh1 * fx * tm3;
        dky = qh1 * fy * tm2;
        dkz = qh1 * fz * tm1;
    end

    kxn = kx + dkx;
    kyn = ky + dky;
    kzn = kz + dkz;

    skx = power(kxn,2);
    sky = power(kyn,2);
    skz = power(kzn,2);

    sk = skx+sky+skz;
    gk = hhm * sk;
    energy = 2*gk/(1+sqrt(1+4*alpha*gk));

end