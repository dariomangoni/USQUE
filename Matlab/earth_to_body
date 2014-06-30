function [vettore_bodyframe] = earth_to_body(vettore_earthframe, quaternione)

quaternione_conj = quaternione.*[1;-1;-1;-1];
vettore_bodyframe = quaternione_moltiplicazione([0;vettore_earthframe],quaternione_conj);
vettore_bodyframe = quaternione_moltiplicazione(quaternione,vettore_bodyframe);
vettore_bodyframe = vettore_bodyframe(2:4);
end
