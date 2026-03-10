function absval = max_abs_diff(C,G,H,delta,Low,Upp,X_inc,spsolver,n)

    Lows = max(Low - X_inc, -delta);
    Upps = min(Upp - X_inc, delta);
    if spsolver == 1 % Stefan's crappy 10line solver
        [~, mdec] = bqmin(-H, -G, Lows, Upps);
        posval = -(-C + mdec);
        [~, mdec] = bqmin(H, G, Lows, Upps);
        negval = -(C + mdec);
    elseif spsolver == 2 % Arnold Neumaier's minq5
        H = (H + H')/2;
        [~, mdec, minq_err] = minqsw(0, -G, -H, Lows', Upps', 0, zeros(n, 1));
        if minq_err < 0
            error('MINQ failed.')
        end
        posval = -(-C + mdec);
        [~, mdec, minq_err] = minqsw(0, G, H, Lows', Upps', 0, zeros(n, 1));
        if minq_err < 0
            error('MINQ failed.')
        end
        negval = -(C + mdec);
    elseif spsolver == 3 % Arnold Neumaier's minq8
        data.gam = 0;
        data.c = -G;
        data.b = zeros(n, 1);
        [tmp1, tmp2] = ldl(-H);
        data.D = diag(tmp2);
        data.A = tmp1';
        [~, mdec] = minq8(data, Lows', Upps', zeros(n, 1), 10 * n);
        posval = -(-C + mdec);
        data.c = G;
        [tmp1, tmp2] = ldl(H);
        data.D = diag(tmp2);
        data.A = tmp1';
        [~, mdec] = minq8(data, Lows', Upps', zeros(n, 1), 10 * n);
        negval = -(C + mdec);
    end

    absval = max(abs(posval),abs(negval));
end