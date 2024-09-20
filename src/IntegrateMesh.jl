function IntegrateMesh(Γ1::CompScienceMeshes.Mesh{3,3,Float64}, Γ2::CompScienceMeshes.Mesh{3,3,Float64}; operator="singlelayer")

    #unpack the mesh
    vertices1 = Γ1.vertices
    faces1 = Γ1.faces
    n1 = length(faces1)

    vertices2 = Γ2.vertices
    faces2 = Γ2.faces
    n2 = length(faces2)

    #iterate over faces.
    A = zeros(n1, n2)
    for i in 1:n1
        for j in 1:n2
            if operator == "singlelayer"
                L0, ~, ~, ~ = GalerkinLaplaceTriGS(
                    vertices2[faces2[j][1]],
                    vertices2[faces2[j][2]],
                    vertices2[faces2[j][3]],
                    vertices1[faces1[i][1]],
                    vertices1[faces1[i][2]],
                    vertices1[faces1[i][3]],
                )
                A[i, j] = L0 / 4 / pi
            elseif operator == "doublelayer"
                ~, M0, ~, ~ = GalerkinLaplaceTriGS(
                    vertices2[faces2[j][1]],
                    vertices2[faces2[j][2]],
                    vertices2[faces2[j][3]],
                    vertices1[faces1[i][1]],
                    vertices1[faces1[i][2]],
                    vertices1[faces1[i][3]],
                )
                A[i, j] = M0 / 4 / pi
            elseif operator == "doublelayer_transposed"
                ~, ~, Ld0, ~ = GalerkinLaplaceTriGS(
                    vertices2[faces2[j][1]],
                    vertices2[faces2[j][2]],
                    vertices2[faces2[j][3]],
                    vertices1[faces1[i][1]],
                    vertices1[faces1[i][2]],
                    vertices1[faces1[i][3]],
                )
                A[i, j] = Ld0 / 4 / pi
            elseif operator == "hypersingular"
                ~, ~, ~, Md0 = GalerkinLaplaceTriGS(
                    vertices2[faces2[j][1]],
                    vertices2[faces2[j][2]],
                    vertices2[faces2[j][3]],
                    vertices1[faces1[i][1]],
                    vertices1[faces1[i][2]],
                    vertices1[faces1[i][3]];
                    L_output=0,
                    Md_output=1,
                )
                A[i, j] = Md0 / 4 / pi
            else
                error("Operator does not exist.")
            end
        end
    end

    return A
end
