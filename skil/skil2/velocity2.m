function v=velocity2(y,L,v0,n)
    v=v0* ( simpson2(0,y,n) / simpson2(0,L,n) );
end
