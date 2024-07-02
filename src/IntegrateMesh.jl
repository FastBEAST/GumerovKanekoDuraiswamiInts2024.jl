function IntegrateMesh(Γ1::CompScienceMeshes.Mesh{3, 3, Float64}, Γ2::CompScienceMeshes.Mesh{3, 3, Float64}; operator = "singlelayer")
    #TODO Add operator input: Helmholtz3D etc. Do this when how to integrate into BEAST is figured out because Helmholtz3D is defined in BEAST.
    #Actually I can go into GL and only calculate what is needed to reduce runtime

    #unpack the mesh
    vertices1 = Γ1.vertices
    faces1 = Γ1.faces
    n1 = length(faces1)

    vertices2 = Γ2.vertices
    faces2 = Γ2.faces
    n2 = length(faces2)

    #iterate over faces.
    # A = Matrix{Float64}(undef,n2,n1)
    # fill!(A,0)
    A = zeros(n1,n2)
    for i = 1:n1
        for j = 1:n2
            # println(vertices[faces[i][1]],vertices[faces[i][2]],vertices[faces[i][3]],vertices[faces[j][1]],vertices[faces[j][2]],vertices[faces[j][3]])
            # if i == 2 && j == 32
            #     error("error is here.")
            # end
            # if testMode == 1
            # println("x:$i, y:$j")
            if operator == "singlelayer"
                L0,~,~,~ = GalerkinLaplaceTriGS(vertices2[faces2[j][1]],vertices2[faces2[j][2]],vertices2[faces2[j][3]],vertices1[faces1[i][1]],vertices1[faces1[i][2]],vertices1[faces1[i][3]])
                A[i,j] = L0/4/pi
            elseif operator == "doublelayer"
                ~,M0,~,~ = GalerkinLaplaceTriGS(vertices2[faces2[j][1]],vertices2[faces2[j][2]],vertices2[faces2[j][3]],vertices1[faces1[i][1]],vertices1[faces1[i][2]],vertices1[faces1[i][3]])
                A[i,j] = M0/4/pi
            elseif operator == "doublelayer normal derivative"
                ~,~,~,Md0 = GalerkinLaplaceTriGS(vertices2[faces2[j][1]],vertices2[faces2[j][2]],vertices2[faces2[j][3]],vertices1[faces1[i][1]],vertices1[faces1[i][2]],vertices1[faces1[i][3]], L_output = 0, Md_output = 1)
                A[i,j] = Md0/4/pi
            else
                error("Operator does not exist.")
            end
        end
    end

    return A
end