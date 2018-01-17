function fi = fi(MOD, VecX, VecY, VecU, ~)
    fi = MOD.fqei_dfqeiXYU([VecX; VecY], 'f', 'i');
end
