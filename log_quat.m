function th = log_quat(quat)
    q12 = quatnormalize(quat);

%             axang = quat2axang(q12);
    th = 2 * acos(q12(:,1));
    th(th > pi) = th(th > pi) - 2*pi;
    u = q12(:,2:4) ./ sqrt(1 - q12(:,1).^2);
    axang = [u, th];

    th = axang(:,1:3) .* axang(:,4);
    th(all(isnan(th), 2),:) = 0;
end