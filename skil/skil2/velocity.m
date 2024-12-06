function v=velocity(y,L,v0,n)
    v=v0* ( simpson(0,y,n) / simpson(0,L,n) );
end
