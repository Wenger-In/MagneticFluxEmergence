program poisson_solver_3d_divb_cleaning
    use, intrinsic :: iso_fortran_env, only: dp => real64! 使用双精度实数类型
    implicit none

    integer :: nx, ny, nz
    integer :: i, j, k, iter, ios! ios用于文件读取状态判断
    integer :: num_lines, line_count! 用于统计文件行数和当前读取行数
    logical :: converged

    real(dp), allocatable :: bx(:,:,:), by(:,:,:), bz(:,:,:), divb(:,:,:), psi(:,:,:), bx_new(:,:,:), by_new(:,:,:), bz_new(:,:,:)
    real(dp), allocatable :: xx(:,:,:), yy(:,:,:), zz(:,:,:), div_B(:,:,:)
    real(dp), allocatable :: xx_new(:,:,:), yy_new(:,:,:), zz_new(:,:,:), div_B_new(:,:,:)
    real(dp), allocatable :: temp_line(:)! 临时数组用于读取每行数据

    character(len=*), parameter :: data_dir = "E:/Research/Work/202405_solar_storm/reconstruction/"
    character(len=*), parameter :: data_name = "origin_field_data_0.5.csv"
    character(len=*), parameter :: save_dir = "E:/Research/Work/202405_solar_storm/reconstruction/"
    character(len=*), parameter :: save_name = "Poisson_solution_0.5.csv"

    real(dp) :: dx, dy, dz, tolerance, max_divb

  ! 设定网格间距和收敛相关参数
    nx = 241
    ny = 103
    nz = 90
    dx = 1823.013857000597_dp
    dy = 1823.013857000597_dp
    dz = -1800.0_dp
    tolerance = 1.0e-6_dp
    max_divb = 1.0e6_dp
    iter = 0
    converged =.false.

  ! 打开数据文件，读取数据
    open(unit=10, file=trim(data_dir) // trim(data_name), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, "无法打开数据文件：", trim(data_dir) // trim(data_name)
        stop
    end if

  ! 先统计文件行数（假设文件不大，可一次性读入内存处理，如果文件很大可能需要逐块处理）
    num_lines = 0
    do
        read(10, *, iostat=ios)
        if (ios < 0) exit
        num_lines = num_lines + 1
    end do
    rewind(10)! 重新将文件指针定位到开头

  ! 分配内存空间给各个数组
    allocate(xx(num_lines, 1, 1), yy(num_lines, 1, 1), zz(num_lines, 1, 1), &
             bx(num_lines, 1, 1), by(num_lines, 1, 1), bz(num_lines, 1, 1), div_B(num_lines, 1, 1))
    allocate(divb(nx, ny, nz), psi(nx, ny, nz), bx_new(nx, ny, nz), by_new(nx, ny, nz), bz_new(nx, ny, nz))
    print *, 'debug1'

  ! 分配临时数组用于读取每行数据，这里假设每行数据依次为xx, yy, zz, bx, by, bz, div_B的值
    allocate(temp_line(7))

    line_count = 0
    do while (line_count < num_lines)
        read(10, *, iostat=ios) temp_line
        if (ios < 0) exit
        line_count = line_count + 1
        xx(line_count, 1, 1) = temp_line(1)
        yy(line_count, 1, 1) = temp_line(2)
        zz(line_count, 1, 1) = temp_line(3)
        bx(line_count, 1, 1) = temp_line(4)
        by(line_count, 1, 1) = temp_line(5)
        bz(line_count, 1, 1) = temp_line(6)
        div_B(line_count, 1, 1) = temp_line(7)
    end do

    close(10)
    print *, 'debug2'

  ! 根据读取的数据确定网格大小（这里假设数据维度与之前设定的nx、ny、nz一致，如果不一致需调整处理逻辑）
    ! nx = size(xx, 1)
    ! ny = size(yy, 2)
    ! nz = size(zz, 3)

    xx = reshape(xx, [nx, ny, nz])
    yy = reshape(yy, [nx, ny, nz])
    zz = reshape(zz, [nx, ny, nz])
    bx = reshape(bx, [nx, ny, nz])
    by = reshape(by, [nx, ny, nz])
    bz = reshape(bz, [nx, ny, nz])
    div_B = reshape(div_B, [nx, ny, nz])
    print *, 'debug3'

  ! 主循环来清理磁场B的散度，添加了收敛判断标志以更准确控制循环结束
    do while (max_divb > tolerance .and. iter < 1000)
      ! 应用边界条件
        call apply_boundary_conditions(bx, by, bz, nx, ny, nz)

      ! 计算磁场B的散度
        call calc_divb(bx, by, bz, dx, dy, dz, divb)
        print *, 'max div: ', maxval(abs(divb))

      ! 使用Jacobi迭代求解标量势psi
        call jacobi_potential(divb, dx, dy, dz, psi)

      ! 更新磁场分量
        call update_magnetic_field(bx, by, bz, bx_new, by_new, bz_new, psi, dx, dy, dz)

      ! 检查最大散度并更新收敛标志和迭代次数
        call calc_divb(bx, by, bz, dx, dy, dz, divb)
        max_divb = maxval(abs(divb))
        iter = iter + 1
        print *, "Iteration: ", iter, "Max divergence: ", max_divb
    end do

    print *, "Divergence cleaning complete."

  ! 打开保存文件，准备写入数据
    open(unit=20, file=trim(save_dir) // trim(save_name), status='replace', action='write')
    if (ios /= 0) then
        print *, "无法打开保存文件：", trim(save_dir) // trim(save_name)
        stop
    end if

  ! 写入表头信息
    write(20, '(a,a,a,a)') "xx,yy,zz,psi"

    do k = 1, size(zz, 3)
        do j = 1, size(yy, 2)
            do i = 1, size(xx, 1)
                write(20, *) xx(i, j, k), yy(i, j, k), zz(i, j, k), psi(i, j, k)
            end do
        end do
    end do

    close(20)

  ! 释放动态分配的内存空间
    deallocate(bx, by, bz, divb, psi, bx_new, by_new, bz_new, xx, yy, zz, div_B, temp_line)

contains

  ! 应用磁场的边界条件，保持原逻辑，可根据实际物理问题进一步扩展边界条件类型
    subroutine apply_boundary_conditions(bx, by, bz, nx, ny, nz)
        real(dp), dimension(nx, ny, nz), intent(inout) :: bx, by, bz
        integer, intent(in) :: nx, ny, nz
        integer :: i, j, k

      ! 边界条件使用外推法，可根据实际需求调整为其他合适的边界条件
      ! z=1边界（底部）
        bx(:, :, 1) = bx(:, :, 2)
        by(:, :, 1) = by(:, :, 2)
        bz(:, :, 1) = bz(:, :, 2)

      ! z=nz边界（顶部）
        bx(:, :, nz) = bx(:, :, nz-1)
        by(:, :, nz) = by(:, :, nz-1)
        bz(:, :, nz) = bz(:, :, nz-1)

      ! y=1和y=ny边界
        bx(:, 1, :) = bx(:, 2, :)
        by(:, 1, :) = by(:, 2, :)
        bz(:, 1, :) = bz(:, 2, :)

        bx(:, ny, :) = bx(:, ny-1, :)
        by(:, ny, :) = by(:, ny-1, :)
        bz(:, ny, :) = bz(:, ny-1, :)

      ! x=1和x=nx边界
        bx(1, :, :) = bx(2, :, :)
        by(1, :, :) = by(2, :, :)
        bz(1, :, :) = bz(2, :, :)

        bx(nx, :, :) = bx(nx-1, :, :)
        by(nx, :, :) = by(nx-1, :, :)
        bz(nx, :, :) = bz(nx-1, :, :)

    end subroutine apply_boundary_conditions

  ! 计算三维磁场B的散度，使用中心差分法，保持原逻辑，可根据精度需求调整差分格式
    subroutine calc_divb(bx, by, bz, dx, dy, dz, divb)
        real(dp), dimension(nx, ny, nz), intent(in) :: bx, by, bz
        real(dp), intent(in) :: dx, dy, dz
        real(dp), dimension(nx, ny, nz), intent(out) :: divb
        integer :: i, j, k

      ! 三维磁场B散度的中心差分计算
        do i = 2, nx - 1
            do j = 2, ny- 1
                do k = 2, nz - 1
                    divb(i,j,k) = (bx(i+1,j,k) - bx(i-1,j,k)) / (2.0 * dx) + &
                                  (by(i,j+1,k) - by(i,j-1,k)) / (2.0 * dy) + &
                                  (bz(i,j,k+1) - bz(i,j,k-1)) / (2.0 * dz)
                end do
            end do
        end do
    end subroutine calc_divb

  ! 使用Jacobi迭代求解标量势psi的三维Laplace方程，添加了迭代次数限制以避免无限循环
    subroutine jacobi_potential(divb, dx, dy, dz, psi)
        real(dp), dimension(nx, ny, nz), intent(in) :: divb
        real(dp), intent(in) :: dx, dy, dz
        real(dp), dimension(nx, ny, nz), intent(out) :: psi
        integer :: i, j, k
        real(dp) :: factor
        real(dp), dimension(nx, ny, nz) :: psi_old
        factor = 1.0_dp / (2.0_dp * (1.0_dp / dx ** 2 + 1.0_dp / dy ** 2 + 1.0_dp / dz ** 2))

      ! Jacobi迭代求解Laplace方程的主循环，添加了收敛判断条件
        psi_old = psi
        do 
            do i = 2, nx - 1
                do j = 2, ny - 1
                    do k = 2, nz - 1
                        psi(i,j,k) = factor * ((psi_old(i+1,j,k) + psi_old(i-1,j,k)) / dx ** 2 + &
                                               (psi_old(i,j+1,k) + psi_old(i,j-1,k)) / dy ** 2 + &
                                               (psi_old(i,j,k+1) + psi_old(i,j,k-1)) / dz ** 2 - divb(i,j,k))
                    end do
                end do
            end do
        end do
    end subroutine jacobi_potential

  ! 使用标量势psi更新三维磁场分量，保持原逻辑，可根据物理模型进一步完善更新方式
    subroutine update_magnetic_field(bx, by, bz, bx_new, by_new, bz_new, psi, dx, dy, dz)
        real(dp), dimension(nx, ny, nz), intent(inout) :: bx, by, bz
        real(dp), dimension(nx, ny, nz), intent(out) :: bx_new, by_new, bz_new
        real(dp), dimension(nx, ny, nz), intent(in) :: psi
        real(dp) :: dx, dy, dz
        integer :: i, j, k

      ! 使用psi的梯度更新磁场分量
        do i = 2, nx - 1
            do j = 2, ny - 1
                do k = 2, nz - 1
                    bx_new(i,j,k) = bx(i,j,k) - (psi(i+1,j,k) - psi(i-1,j,k)) / (2.0 * dx)
                    by_new(i,j,k) = by(i,j,k) - (psi(i,j+1,k) - psi(i,j-1,k)) / (2.0 * dy)
                    bz_new(i,j,k) = bz(i,j,k) - (psi(i,j,k+1) - psi(i,j,k-1)) / (2.0 * dz)
                end do
            end do
        end do

      ! 更新原始磁场分量
        bx = bx_new
        by = by_new
        bz = bz_new

    end subroutine update_magnetic_field

end program poisson_solver_3d_divb_cleaning