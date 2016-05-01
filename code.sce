// -------
// Question 15
// -------

lambda01 = 2.40483
lambda11 = 3.83171
lambda03 = 8.65373

// t = 0
function z = init_positions()
  z = zeros(N_eta, N_theta + 1)

  for i = 1:N_eta
    z(i,:) = z(i,:) + besselj(0, lambda03 * i * d_eta)
  end
endfunction

// i > 1
function c = get_coeff(mat, i, j)
  // Conditions aux limites
  if (i > N_eta) then
    c = 0
    return
  end

  // Conditions periodiques en theta
  if (j < 1) then
    j = (N_theta + 1) + j
  elseif (j > N_theta + 1) then
    j = modulo(j, N_theta + 1)
  end

  c = mat(i, j)
endfunction

function nxt = next_positions(cur, prev)
  nxt = zeros(N_eta, N_theta + 1)

  // eta = 0, i.e. i = 1
  s = 0
  for j = 1:(N_theta+1)
    s = s + cur(2, j) - cur(1, 1)
  end
  nxt(1,:) = 2 * cur(1, 1) - prev(1, 1) + d_tau^2 * 4/(d_eta)^2 * 1/(N_theta + 1) * s

  // eta > 0, i.e. i > 1
  for i = 2:N_eta
    for j = 1:(N_theta + 1)
      nxt(i, j) =  1/(i * d_eta) * (get_coeff(cur, i+1, j) - get_coeff(cur, i-1, j))/(2 * d_eta)
      nxt(i, j) = nxt(i, j) + (get_coeff(cur, i+1, j) - 2 * get_coeff(cur, i, j) + get_coeff(cur, i-1, j))/(d_eta)^2
      nxt(i, j) = nxt(i, j) + 1/(i * d_eta)^2 * (get_coeff(cur, i, j+1) - 2 * get_coeff(cur, i, j) + get_coeff(cur, i, j-1))/(d_theta)^2
      nxt(i, j) = nxt(i, j) * d_tau^2
      nxt(i, j) = nxt(i, j) + 2 * get_coeff(cur, i, j) - get_coeff(prev, i, j)
    end
  end
endfunction

function plot_positions(z)
  eta = linspace(0, 1, N_eta)
  theta = linspace(0, 1, N_theta + 1)
  plot3d(eta, theta, z)
endfunction

function q15_animation()
  cur = init_positions()
  prev = cur // Conditions aux limites sur la derivee
  plot_positions(cur)

  N_tau = 100
  d_tau = 0.01
  N_theta = 80
  d_theta = 1 / (N_theta + 1) // La matrice est de taille N_eta * (N_theta + 1)
  N_eta = 40
  d_eta = 1 / N_eta
  c = 1

  for t = 1:N_tau
    drawlater;
    nxt = next_positions(cur, prev)
    clf();
    ax = gca()
    ax.data_bounds = [0, 0, -1 ; 1, 1, 1]
    xtitle('Numérique', 'eta', 'theta', 'w')
    plot_positions(nxt)
    drawnow;
    prev = cur
    cur = nxt
  end
endfunction
