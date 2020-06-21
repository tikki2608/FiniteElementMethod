function [mesh] = make_project2_mesh(h, r0, r1, etype)
    nex = ceil((r1 - r0) / h);
    ney = ceil((pi/2 * r1) / h);
    
    if strcmp(etype, 'q4')
        mesh = mesh_unit_square_q4(nex, ney);
    elseif strcmp(etype, 't3')
        mesh = mesh_unit_square_q4(nex, ney);
        mesh = convert_q4_to_t3(mesh);
    elseif strcmp(etype, 'q8')
        mesh = mesh_unit_square_q8(nex, ney);
    elseif strcmp(etype, 'q9')
        mesh = mesh_unit_square_q9(nex, ney);
    elseif strcmp(etype, 't6')
        mesh = mesh_unit_square_q9(nex, ney);
        mesh = convert_q9_to_t6(mesh);
    end
    mesh = transform_to_ring(mesh, r0, r1);
    fprintf('Generated mesh with %d elements and %d nodes.\n', ...
            size(mesh.conn,2), size(mesh.x,2));
end

function [mesh] = mesh_unit_square_q9(nex, ney)
    mesh.x = zeros(2, (nex+1)*ney+1);
    mesh.conn = zeros(9, nex*ney);
    nnx = 2*nex + 1;
    nny = 2*ney + 1;   
    for j=1:nny
        for i=1:nnx
            mesh.x(:, i+(j-1)*nnx) = [(i-1)/(nnx-1); (j-1)/(nny-1)];            
        end
    end
    for j=1:ney
        for i=1:nex
            n1 = 2*i-1 + (j-1)*nnx*2;
            n4 = n1+2*nnx;
            n8 = n1 + nnx;
            mesh.conn(:,i+(j-1)*nex) = [n1; n1+2; n4+2; n4;
                                        n1+1; n8+2; n4+1; n8; n8 + 1];              
        end
    end
    % Splits each Q9 element into 4 Q4 elements (just for plotting).
    mesh.pconn = [mesh.conn([1,5,9,8], :), ...
                  mesh.conn([5,2,6,9], :), ... 
                  mesh.conn([9,6,3,7], :), ... 
                  mesh.conn([8,9,7,4], :)];
end


function [mesh] = convert_q9_to_t6(mesh)
    conn = mesh.conn;
    % Splits each Q9 element into 2 T6 elements.
    mesh.conn = [conn([1;2;3;5;6;9],:), conn([3;4;1;7;8;9],:)];
    % Splits each T6 element into 4 T3 elements (just for plotting).
    mesh.pconn = [mesh.conn([1,4,6],:),...
                  mesh.conn([2,6,4],:),...
                  mesh.conn([2,6,5],:),...
                  mesh.conn([5,6,3],:)];
end


function [mesh] = mesh_unit_square_q8(nex, ney)
    % Number of nodes along odd rows.    
    nx1 = 1 + 2*nex;
    % Number of nodes along even rows.
    nx2 = 1 + nex;
    % Total number of nodes.
    nn = (1+ney) * nx1 + ney*nx2;    
    mesh.x = zeros(2, nn);
    mesh.conn = zeros(8, nex*ney);
    % Defines bottom row of nodes.
    mesh.x(1, 1:nx1) = linspace(0, 1, nx1);
    mesh.x(2, 1:nx1) = 0.0;
    ct = nx1;
    for j=1:ney
        J = ct + 1;
        % For each row of elements define middle row of nodes.
        mesh.x(1, ct+1:ct+nx2) = linspace(0, 1, nx2);
        mesh.x(2, ct+1:ct+nx2) = j/ney - 0.5/ney;
        ct = ct + nx2;
        % For each row of elements define top row of nodes.
        mesh.x(1, ct+1:ct+nx1) = linspace(0, 1, nx1);
        mesh.x(2, ct+1:ct+nx1) = j/ney;
        ct = ct + nx1;        
        for i = 1:nex
            I = 2*i-1 + (j-1)*(nx1 + nx2);
            e = i+nex*(j-1);
            mesh.conn(:, e) = [I; I+2; I+nx1+nx2+2; I+nx1+nx2;
                               I+1; J+i; I+nx1+nx2+1; J+i-1];
        end        
    end
    mesh.pconn = mesh.conn([1;5;2;6;3;7;4;8], :);
end

function [mesh] = mesh_unit_square_q4(nex, ney)
    mesh.x = zeros(2, (nex+1)*ney+1);
    mesh.conn = zeros(4, nex*ney);
    for j=1:ney+1
        for i=1:nex+1
            I = i+(j-1)*(nex+1);
            mesh.x(:, I) = [(i-1)/nex; (j-1)/ney];
            if i <= nex && j <= ney
                mesh.conn(:,i+(j-1)*nex) = [I; I+1; I+1+nex+1; I+nex+1];
            end
        end
    end
end

function [mesh] = convert_q4_to_t3(mesh)
    conn = mesh.conn;
    mesh.conn = [conn([1;2;3],:), conn([3;4;1],:)];
end

function [mesh] = transform_to_ring(mesh, r0, r1)
    % Map from square to quarter of the ring.
    r = r0 + mesh.x(1,:)*(r1 - r0);
    q = mesh.x(2,:) * pi / 2;    
    mesh.x(1,:) = cos(q) .* r;
    mesh.x(2,:) = sin(q) .* r;    
end