driver("Rec")

lambda01 = 2.40483
lambda11 = 3.83171
lambda03 = 8.65373

// -------
// Question 15
// -------

function z = init_positions1()
  z = zeros(N_eta, N_theta + 1)

  for i = 1:N_eta
    z(i,:) = besselj(0, lambda03 * i * d_eta)
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

function plot_positions(z, plot_title)
  eta = linspace(0, 1, N_eta)
  theta = linspace(0, 2 * %pi, N_theta + 1)
  
  x = zeros(N_eta, N_theta + 1)
  y = zeros(N_eta, N_theta + 1)

  for i = 1:N_eta
    for j = 1:(N_theta + 1)
      x(i, j) = eta(i) * cos(theta(j))
      y(i, j) = eta(i) * sin(theta(j))
    end
  end
  
  ax = gca()
  ax.data_bounds = [-1, -1, -1 ; 1, 1, 1]
  xtitle(plot_title, 'x', 'y', 'w')
  surf(x, y, z)
endfunction

function num_animation(init_positions)
  cur = init_positions()
  prev = cur // Conditions aux limites sur la derivee
  plot_positions(cur, 'Numérique')

  for t = 1:N_tau
    drawlater;
    nxt = next_positions(cur, prev)
    clf();
    plot_positions(nxt, 'Numérique')
    drawnow;
    im_name = 'image' + string(t) + '.gif'
    xs2gif(0, im_name)
    prev = cur
    cur = nxt
  end
endfunction

function q15_animation()
  N_tau = 100
  d_tau = 0.01
  N_theta = 80
  d_theta = 1 / (N_theta + 1) // La matrice est de taille N_eta * (N_theta + 1)
  N_eta = 40
  d_eta = 1 / N_eta

  num_animation(init_positions1)
endfunction

// -------
// Question 16
// -------

function z = get_positions_ex(n)
  z = zeros(N_eta, N_theta + 1)

  for i = 1:N_eta
    z(i,:) = cos(lambda03 * n * d_tau) * besselj(0, lambda03 * i * d_eta)
  end
endfunction

function q16_animation()
  N_tau = 100
  d_tau = 0.01
  N_theta = 80
  d_theta = 1 / (N_theta + 1) // La matrice est de taille N_eta * (N_theta + 1)
  N_eta = 40
  d_eta = 1 / N_eta

  cur = init_positions1()
  prev = cur // Conditions aux limites sur la derivee

  for t = 1:N_tau
    drawlater;
    nxt = next_positions(cur, prev)
    clf();
    subplot(1, 2, 1)
    plot_positions(nxt, 'Numérique')
    subplot(1, 2, 2)
    plot_positions(get_positions_ex(t), 'Exacte')
    drawnow;
    im_name = 'image' + string(t) + '.gif'
    xs2gif(0, im_name)
    prev = cur
    cur = nxt
  end
endfunction

function relerror(ex_positions, init_positions, ttl)
  CFL_start = 0.2
  CFL_step = 0.2

  for CFL = CFL_start:CFL_step:1
    cur = init_positions()
    prev = cur

    d_tau = CFL * d_eta * d_theta
    err = zeros(1, N_tau)

    for t = 1:N_tau
      nxt = next_positions(cur, prev)
      prev = cur
      cur = nxt

      num = nxt(1, 1)
      ex = ex_positions(t)
      ex = ex(1, 1)

      err(1, t) = (num - ex)/abs(ex)
    end

    subplot(2, 3, CFL / CFL_step)
    xtitle('CFL = ' + string(CFL), 'Temps', 'Erreur relative')
    plot2d(1:N_tau, err)
  end
endfunction

function q16_relerror()
  N_tau = 100
  N_theta = 80
  d_theta = 1 / (N_theta + 1) // La matrice est de taille N_eta * (N_theta + 1)
  N_eta = 40
  d_eta = 1 / N_eta

  relerror(get_positions_ex, init_positions1, 'q16relerror')
endfunction

// -------
// Question 17
// -------

function z = get_positions_ex2(n)
  z = zeros(N_eta, N_theta + 1)

  for i = 1:N_eta
    for j = 1:N_theta
      z(i, j) = cos(lambda01 * n * d_tau) * besselj(0, lambda01 * i * d_eta)
      z(i, j) = z(i, j) + cos(lambda11 * n * d_tau) * cos(j * d_theta) * besselj(1, lambda11 * i * d_eta)
    end
  end
endfunction

function z = init_positions2()
  z = zeros(N_eta, N_theta + 1)

  for i = 1:N_eta
    for j = 1:(N_theta + 1)
      z(i, j) = besselj(0, lambda01 * i * d_eta) * besselj(1, lambda11 * i * d_eta) * cos(j * d_theta)
    end
  end
endfunction

function q17_relerror()
  N_tau = 100
  N_theta = 80
  d_theta = 1 / (N_theta + 1) // La matrice est de taille N_eta * (N_theta + 1)
  N_eta = 40
  d_eta = 1 / N_eta

  relerror(get_positions_ex2, init_positions2, 'q17relerror')
endfunction

// -------
// Question 19
// -------

function q19_relerror()
  CFL = 0.8
  N_tau = 100
  grids = [80, 40 ; 40, 20 ; 160, 80]

  errs = zeros(1, N_tau)

  for i = 1:3
    grid = grids(i,:)
    N_theta = grid(1)
    d_theta = 1 / (N_theta + 1)
    N_eta = grid(2)
    d_eta = 1 / N_eta
    d_tau = CFL * d_eta * d_theta

    cur = init_positions1()
    prev = cur

    for t = 1:N_tau
      nxt = next_positions(cur, prev)
      ex = get_positions_ex(t)

      errs(t) = max(abs(nxt - ex))

      prev = cur
      cur = nxt
    end

    xtitle('', 'Temps', 'Erreur relative globale')
    plot2d(1:N_tau, errs, i)
  end

  legend(['(80, 40)';'(40, 20)';'(160, 80)'])
  xs2pdf(0, 'q19_relerror')
endfunction

// -------
// Question 20
// -------

function z = init_positions3()
  z = zeros(N_eta, N_theta + 1)

  for i = 1:N_eta
    for j = 1:(N_theta + 1)
      z(i, j) = besselj(1, lambda11 * i * d_eta) * cos(j * d_theta) / 2
    end
  end
endfunction

function q20_animation()
  N_tau = 200
  d_tau = 0.01
  N_theta = 80
  d_theta = 1 / (N_theta + 1) // La matrice est de taille N_eta * (N_theta + 1)
  N_eta = 40
  d_eta = 1 / N_eta

  num_animation(init_positions3)
endfunction
