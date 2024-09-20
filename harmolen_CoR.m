function C=harmolen_CoR(v_in,p,angle)
    % coefficient of restitution 
   
    if nargin>2
        % this is a crest hit
        C = [(1-p.k_txc*abs(v_in(1,:)).^(3/5)).^.5; (1-p.k_tyc*abs(v_in(2,:)).^(3/5)).^.5; (1-p.k_nc*abs(v_in(3,:)).^(3/5)).^.5]; 
    else
        C = [repmat((1-p.k_tf*sqrt(v_in(1,:).^2+v_in(2,:).^2).^(3/5)).^.5,2,1); (1-p.k_nf*abs(v_in(3,:)).^(3/5)).^.5]; 
    end
end