! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module sw_test_cases

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_log 


   contains


   subroutine setup_sw_test_case(domain)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Configure grid metadata and model state for the shallow water test case 
   !   specified in the namelist
   !
   ! Output: block - a subset (not necessarily proper) of the model domain to be
   !                 initialized
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (domain_type), intent(inout) :: domain

      integer :: i
      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: statePool, meshPool

      integer, pointer :: config_test_case

      call mpas_pool_get_config(domain % configs, 'config_test_case', config_test_case)

      if (config_test_case == 0) then
         call mpas_log_write('Using initial conditions supplied in input file')

         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

            ! Compute xEdgeMid, yEdgeMid, zEdgeMid, xEdgeMidSphere, yEdgeMidSphere and zEdgeMidSphere using MPAS
            call xyzEdgeMid(meshPool)

            ! Compute latEdgeMid, lonEdgeMid using MPAS
            call computeLatLonEdgeMid(meshPool)

            ! Compute edgeSignOnCell and edgeSignOnVertex 
            call ocn_init_routines_setup_sign_and_index_fields(meshPool)

            ! Compute weightsOnEdge in the baricentric way using MPAS (not the external C code)
            call compute_weightsOnEdge_baricentric(meshPool)

            ! Compute the weights for h_vertex and pv_vertex needed in the HCm discretization 
            call weightsPV(meshPool)

            ! Compute the weights for h_edge needed in the HCm discretization 
            call weightsHEDGE(meshPool)

            call sw_test_case_0(meshPool, statePool)

            call mpas_pool_initialize_time_levels(statePool)

            block_ptr => block_ptr % next
         end do

      else if (config_test_case == 1) then
         call mpas_log_write('Setting up shallow water test case 1')
         call mpas_log_write(' -- Advection of Cosine Bell over the Pole')

         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

            call sw_test_case_1(meshPool, statePool)

            call mpas_pool_initialize_time_levels(statePool)

            block_ptr => block_ptr % next
         end do

      else if (config_test_case == 2) then
         call mpas_log_write('Setting up shallow water test case 2')
         call mpas_log_write(' -- Setup shallow water test case 2: Global Steady State Nonlinear Zonal Geostrophic Flow')

         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

            ! Compute xEdgeMid, yEdgeMid, zEdgeMid, xEdgeMidSphere, yEdgeMidSphere and zEdgeMidSphere using MPAS
            call xyzEdgeMid(meshPool)
         
            ! Compute latEdgeMid, lonEdgeMid using MPAS
            call computeLatLonEdgeMid(meshPool)
            
            ! Compute edgeSignOnCell and edgeSignOnVertex 
            call ocn_init_routines_setup_sign_and_index_fields(meshPool) 
            
            ! Compute weightsOnEdge in the baricentric way using MPAS (not the external C code)
            call compute_weightsOnEdge_baricentric(meshPool)

            ! Compute the weights for h_vertex and pv_vertex needed in the HCm discretization 
            call weightsPV(meshPool)

            ! Compute the weights for h_edge needed in the HCm discretization 
            call weightsHEDGE(meshPool)

            call sw_test_case_2(meshPool, statePool)

            call mpas_pool_initialize_time_levels(statePool)

            block_ptr => block_ptr % next
         end do

      else if (config_test_case == 5) then
         call mpas_log_write('Setting up shallow water test case 5')
         call mpas_log_write(' -- Setup shallow water test case 5: Zonal Flow over an Isolated Mountain')

         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

            ! Compute xEdgeMid, yEdgeMid, zEdgeMid, xEdgeMidSphere,
            ! yEdgeMidSphere and zEdgeMidSphere using MPAS
            call xyzEdgeMid(meshPool)

            ! Compute latEdgeMid, lonEdgeMid using MPAS
            call computeLatLonEdgeMid(meshPool)

            ! Compute edgeSignOnCell and edgeSignOnVertex 
            call ocn_init_routines_setup_sign_and_index_fields(meshPool)

            ! Compute weightsOnEdge in the baricentric way using MPAS (not the external C code)
            call compute_weightsOnEdge_baricentric(meshPool)

            ! Compute the weights for h_vertex and pv_vertex needed in the HCm discretization 
            call weightsPV(meshPool)

            ! Compute the weights for h_edge needed in the HCm discretization 
            call weightsHEDGE(meshPool)

            call sw_test_case_5(meshPool, statePool)

            call mpas_pool_initialize_time_levels(statePool)

            block_ptr => block_ptr % next
         end do

      else if (config_test_case == 6) then
         call mpas_log_write('Setting up shallow water test case 6')
         call mpas_log_write(' -- Rossby-Haurwitz Wave')

         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

            ! Compute xEdgeMid, yEdgeMid, zEdgeMid, xEdgeMidSphere,
            ! yEdgeMidSphere and zEdgeMidSphere using MPAS
            call xyzEdgeMid(meshPool)

            ! Compute latEdgeMid, lonEdgeMid using MPAS
            call computeLatLonEdgeMid(meshPool)

            ! Compute edgeSignOnCell and edgeSignOnVertex 
            call ocn_init_routines_setup_sign_and_index_fields(meshPool)

            ! Compute weightsOnEdge in the baricentric way using MPAS (not the external C code)
            call compute_weightsOnEdge_baricentric(meshPool)

            ! Compute the weights for h_vertex and pv_vertex needed in the HCm discretization 
            call weightsPV(meshPool)

            ! Compute the weights for h_edge needed in the HCm discretization 
            call weightsHEDGE(meshPool)

            call sw_test_case_6(meshPool, statePool)

            call mpas_pool_initialize_time_levels(statePool)

            block_ptr => block_ptr % next
         end do

      else if (config_test_case == 7) then  
         !test case on planar mesh, created to test the new operators for HCm
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

            ! Compute xEdgeMid, yEdgeMid, zEdgeMid, xEdgeMidSphere, yEdgeMidSphere and zEdgeMidSphere using MPAS
            call xyzEdgeMid(meshPool)

            call sw_test_case_7(meshPool, statePool)

            call mpas_pool_initialize_time_levels(statePool)

            block_ptr => block_ptr % next
         end do

      else
         call mpas_log_write('Only test case 1, 2, 5, 6 and 7 (a planar case) are currently supported.')
         stop
      end if

   end subroutine setup_sw_test_case


   subroutine sw_test_case_0(meshPool, statePool)
   
      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool
      type (mpas_pool_type), intent(inout) :: statePool

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex, xEdge, yEdge, zEdge, dvEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: xEdgeMid, yEdgeMid, zEdgeMid, xEdgeMidSphere, yEdgeMidSphere, zEdgeMidSphere
      real (kind=RKIND), dimension(:), pointer :: areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex, areaPV, areaHEdge

      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
      call mpas_pool_get_array(meshPool, 'xEdgeMid', xEdgeMid)
      call mpas_pool_get_array(meshPool, 'yEdgeMid', yEdgeMid)
      call mpas_pool_get_array(meshPool, 'zEdgeMid', zEdgeMid)
      call mpas_pool_get_array(meshPool, 'xEdgeMidSphere', xEdgeMidSphere)
      call mpas_pool_get_array(meshPool, 'yEdgeMidSphere', yEdgeMidSphere)
      call mpas_pool_get_array(meshPool, 'zEdgeMidSphere', zEdgeMidSphere)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(meshPool, 'areaPV', areaPV)
      call mpas_pool_get_array(meshPool, 'areaHEdge', areaHEdge)

      !xCell = xCell * a
      !yCell = yCell * a
      !zCell = zCell * a
      !xVertex = xVertex * a
      !yVertex = yVertex * a
      !zVertex = zVertex * a
      !xEdge = xEdge * a
      !yEdge = yEdge * a
      !zEdge = zEdge * a
      !xEdgeMid = xEdgeMid * a
      !yEdgeMid = yEdgeMid * a
      !zEdgeMid = zEdgeMid * a
      !xEdgeMidSphere = xEdgeMidSphere * a
      !yEdgeMidSphere = yEdgeMidSphere * a
      !zEdgeMidSphere = zEdgeMidSphere * a
      !dvEdge = dvEdge * a
      !dcEdge = dcEdge * a
      !areaCell = areaCell * a**2.0
      !areaTriangle = areaTriangle * a**2.0
      !kiteAreasOnVertex = kiteAreasOnVertex * a**2.0
      !areaPV = areaPV * a**2.0
      !areaHEdge = areaHEdge * a**2.0

   end subroutine sw_test_case_0


   subroutine sw_test_case_1(meshPool, statePool)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup shallow water test case 1: Advection of Cosine Bell over the Pole
   !
   ! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical 
   !            Approximations to the Shallow Water Equations in Spherical 
   !            Geometry" J. of Comp. Phys., 102, pp. 211--224
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool
      type (mpas_pool_type), intent(inout) :: statePool

      real (kind=RKIND), parameter :: u0 = 2.0 * pii * a / (12.0 * 86400.0)
      real (kind=RKIND), parameter :: h0 = 1000.0
      real (kind=RKIND), parameter :: theta_c = 0.0
      real (kind=RKIND), parameter :: lambda_c = 3.0*pii/2.0
      real (kind=RKIND), parameter :: alpha = pii/4.0

      integer :: iCell, iEdge, iVtx
      integer, pointer :: nVertices, nEdges, nCells
      real (kind=RKIND) :: r, v
      real (kind=RKIND), allocatable, dimension(:) :: psiVertex

      integer, dimension(:,:), pointer :: verticesOnEdge

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex, xEdge, yEdge, zEdge, dvEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex
      real (kind=RKIND), dimension(:), pointer :: areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex

      real (kind=RKIND), dimension(:,:), pointer :: u, h

      !
      ! Scale all distances and areas from a unit sphere to one with radius a
      !

      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
      call mpas_pool_get_array(meshPool, 'lonVertex', lonVertex)

      call mpas_pool_get_array(statePool, 'u', u, 1)
      call mpas_pool_get_array(statePool, 'h', h, 1)

      xCell = xCell * a
      yCell = yCell * a
      zCell = zCell * a
      xVertex = xVertex * a
      yVertex = yVertex * a
      zVertex = zVertex * a
      xEdge = xEdge * a
      yEdge = yEdge * a
      zEdge = zEdge * a
      dvEdge = dvEdge * a
      dcEdge = dcEdge * a
      areaCell = areaCell * a**2.0
      areaTriangle = areaTriangle * a**2.0
      kiteAreasOnVertex = kiteAreasOnVertex * a**2.0

      !
      ! Initialize wind field
      !
      allocate(psiVertex(nVertices))
      do iVtx = 1, nVertices
         psiVertex(iVtx) = -a * u0 * ( &
                                       sin(latVertex(iVtx)) * cos(alpha) - &
                                       cos(lonVertex(iVtx)) * cos(latVertex(iVtx)) * sin(alpha) &
                                     )
      end do
      do iEdge = 1, nEdges
         u(1,iEdge) = -1.0 * ( &
                                               psiVertex(verticesOnEdge(2,iEdge)) - &
                                               psiVertex(verticesOnEdge(1,iEdge)) &
                                             ) / dvEdge(iEdge)
      end do
      deallocate(psiVertex)

      !
      ! Initialize cosine bell at (theta_c, lambda_c)
      !
      do iCell = 1, nCells
         r = sphere_distance(theta_c, lambda_c, latCell(iCell), lonCell(iCell), a) 
         if (r < a/3.0) then
            h(1,iCell) = (h0 / 2.0) * (1.0 + cos(pii*r*3.0/a))
         else
            h(1,iCell) = h0 / 2.0
         end if
      end do

   end subroutine sw_test_case_1


   subroutine sw_test_case_2(meshPool, statePool)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup shallow water test case 2: Global Steady State Nonlinear Zonal 
   !                                  Geostrophic Flow
   !
   ! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical 
   !            Approximations to the Shallow Water Equations in Spherical 
   !            Geometry" J. of Comp. Phys., 102, pp. 211--224
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool
      type (mpas_pool_type), intent(inout) :: statePool

      real (kind=RKIND), parameter :: u0 = 2.0 * pii * a / (12.0 * 86400.0)
      real (kind=RKIND), parameter :: gh0 = 29400.0
      real (kind=RKIND), parameter :: alpha = 0.0

      integer :: iCell, iEdge, iVtx, i, k, HCm 
      integer, pointer :: nVertices, nEdges, nCells, nVertLevels
      real (kind=RKIND) :: v
      real (kind=RKIND), allocatable, dimension(:) :: psiVertex

      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: verticesOnEdge, edgesOnCell

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex, xEdge, yEdge, zEdge, dvEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: xEdgeMid, yEdgeMid, zEdgeMid, xEdgeMidSphere, yEdgeMidSphere, zEdgeMidSphere
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex
      real (kind=RKIND), dimension(:), pointer :: latEdge, lonEdge
      real (kind=RKIND), dimension(:), pointer :: latEdgeMid, lonEdgeMid
      real (kind=RKIND), dimension(:), pointer :: areaCell, areaTriangle, fEdge, fVertex 
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex, areaPV, areaHEdge, keInitial

      real (kind=RKIND), dimension(:,:), pointer :: u, h


      !
      ! Scale all distances and areas from a unit sphere to one with radius a
      !

      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
      call mpas_pool_get_array(meshPool, 'xEdgeMid', xEdgeMid)
      call mpas_pool_get_array(meshPool, 'yEdgeMid', yEdgeMid)
      call mpas_pool_get_array(meshPool, 'zEdgeMid', zEdgeMid)
      call mpas_pool_get_array(meshPool, 'xEdgeMidSphere', xEdgeMidSphere)
      call mpas_pool_get_array(meshPool, 'yEdgeMidSphere', yEdgeMidSphere)
      call mpas_pool_get_array(meshPool, 'zEdgeMidSphere', zEdgeMidSphere)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
      call mpas_pool_get_array(meshPool, 'fVertex', fVertex)
      call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(meshPool, 'areaPV', areaPV)
      call mpas_pool_get_array(meshPool, 'areaHEdge', areaHEdge)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
      call mpas_pool_get_array(meshPool, 'lonVertex', lonVertex)
      call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
      call mpas_pool_get_array(meshPool, 'lonEdge', lonEdge) 
      call mpas_pool_get_array(meshPool, 'latEdgeMid', latEdgeMid)
      call mpas_pool_get_array(meshPool, 'lonEdgeMid', lonEdgeMid)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)

      call mpas_pool_get_array(statePool, 'u', u, 1)
      call mpas_pool_get_array(statePool, 'h', h, 1)
      call mpas_pool_get_array(statePool, 'keInitial', keInitial)


      xCell = xCell * a
      yCell = yCell * a
      zCell = zCell * a
      xVertex = xVertex * a
      yVertex = yVertex * a
      zVertex = zVertex * a
      xEdge = xEdge * a
      yEdge = yEdge * a
      zEdge = zEdge * a
      xEdgeMid = xEdgeMid * a
      yEdgeMid = yEdgeMid * a
      zEdgeMid = zEdgeMid * a
      xEdgeMidSphere = xEdgeMidSphere * a
      yEdgeMidSphere = yEdgeMidSphere * a
      zEdgeMidSphere = zEdgeMidSphere * a      
      dvEdge = dvEdge * a
      dcEdge = dcEdge * a
      areaCell = areaCell * a**2.0
      areaTriangle = areaTriangle * a**2.0
      kiteAreasOnVertex = kiteAreasOnVertex * a**2.0
      areaPV = areaPV * a**2.0
      areaHEdge = areaHEdge * a**2.0
      
      !HCm = 0.0 !TRiSK
      HCm = 1.0 

      !
      ! Initialize wind field
      !
      allocate(psiVertex(nVertices))
      do iVtx = 1, nVertices
         psiVertex(iVtx) = -a * u0 * ( &
                                 sin(latVertex(iVtx)) * cos(alpha) - &
                                 cos(lonVertex(iVtx)) * cos(latVertex(iVtx)) * sin(alpha) &
                               )
      end do
      do iEdge = 1,nEdges
         u(1,iEdge) = -1.0 * ( &
                                 psiVertex(verticesOnEdge(2,iEdge)) - &
                                 psiVertex(verticesOnEdge(1,iEdge)) &
                               ) / dvEdge(iEdge)
      end do
      deallocate(psiVertex)

      !
      ! Generate rotated Coriolis field
      !
      if (HCm == 0) then 
         do iEdge = 1, nEdges
            fEdge(iEdge) = 2.0 * omega * &
                               ( -cos(lonEdge(iEdge)) * cos(latEdge(iEdge)) * sin(alpha) + &
                                 sin(latEdge(iEdge)) * cos(alpha) &
                               )
         end do 
      elseif (HCm == 1) then
         do iEdge = 1, nEdges     
            fEdge(iEdge) = 2.0 * omega * &
                               ( -cos(lonEdgeMid(iEdge)) * cos(latEdgeMid(iEdge)) * sin(alpha) + &
                                 sin(latEdgeMid(iEdge)) * cos(alpha) &
                               )
         end do
      end if 

      do iVtx = 1, nVertices
         fVertex(iVtx) = 2.0 * omega * &
                                         (-cos(lonVertex(iVtx)) * cos(latVertex(iVtx)) * sin(alpha) + &
                                          sin(latVertex(iVtx)) * cos(alpha) &
                                         )
      end do

      !
      ! Initialize height field (actually, fluid thickness field)
      !
      do iCell = 1, nCells
         h(1,iCell) = (gh0 - (a * omega * u0 + 0.5 * u0**2.0) * &
                                             (-cos(lonCell(iCell)) * cos(latCell(iCell)) * sin(alpha) + &
                                              sin(latCell(iCell)) * cos(alpha) &
                                             )**2.0 &
                                      ) / &
                                      gravity
      end do

      !
      ! Compute kinetic energy at the initial time in each cell
      !
      keInitial(:,:) = 0.0
      do iCell = 1, nCells
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)
            do k = 1, nVertLevels
               keInitial(k,iCell) = keInitial(k,iCell) + 0.25 * dcEdge(iEdge) * dvEdge(iEdge) * u(k,iEdge)**2.0
            end do
         end do
         do k = 1, nVertLevels
            keInitial(k,iCell) = keInitial(k,iCell) / areaCell(iCell)
         end do
      end do

   end subroutine sw_test_case_2


   subroutine sw_test_case_5(meshPool, statePool)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup shallow water test case 5: Zonal Flow over an Isolated Mountain
   !
   ! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical 
   !            Approximations to the Shallow Water Equations in Spherical 
   !            Geometry" J. of Comp. Phys., 102, pp. 211--224
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool
      type (mpas_pool_type), intent(inout) :: statePool

      real (kind=RKIND), parameter :: u0 = 20.
      real (kind=RKIND), parameter :: gh0 = 5960.0*gravity
      real (kind=RKIND), parameter :: hs0 = 2000.
      real (kind=RKIND), parameter :: theta_c = pii/6.0
      real (kind=RKIND), parameter :: lambda_c = 3.0*pii/2.0
      real (kind=RKIND), parameter :: rr = pii/9.0
      real (kind=RKIND), parameter :: alpha = 0.0

      integer :: iCell, iEdge, iVtx, nTracers, HCm 
      integer, pointer :: nVertices, nEdges, nCells
      real (kind=RKIND) :: r, v
      real (kind=RKIND), allocatable, dimension(:) :: psiVertex

      integer, dimension(:,:), pointer :: verticesOnEdge

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex, xEdge, yEdge, zEdge, dvEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: xEdgeMid, yEdgeMid, zEdgeMid, xEdgeMidSphere, yEdgeMidSphere, zEdgeMidSphere
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex
      real (kind=RKIND), dimension(:), pointer :: latEdge, lonEdge
      real (kind=RKIND), dimension(:), pointer :: latEdgeMid, lonEdgeMid
      real (kind=RKIND), dimension(:), pointer :: areaCell, areaTriangle, h_s, fEdge, fVertex
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex, areaPV, areaHEdge

      real (kind=RKIND), dimension(:,:), pointer :: u, h
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers

      !
      ! Scale all distances and areas from a unit sphere to one with radius a
      !

      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
      call mpas_pool_get_array(meshPool, 'xEdgeMid', xEdgeMid)
      call mpas_pool_get_array(meshPool, 'yEdgeMid', yEdgeMid)
      call mpas_pool_get_array(meshPool, 'zEdgeMid', zEdgeMid)
      call mpas_pool_get_array(meshPool, 'xEdgeMidSphere', xEdgeMidSphere)
      call mpas_pool_get_array(meshPool, 'yEdgeMidSphere', yEdgeMidSphere)
      call mpas_pool_get_array(meshPool, 'zEdgeMidSphere', zEdgeMidSphere)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(meshPool, 'areaPV', areaPV)
      call mpas_pool_get_array(meshPool, 'areaHEdge', areaHEdge)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
      call mpas_pool_get_array(meshPool, 'lonVertex', lonVertex)
      call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
      call mpas_pool_get_array(meshPool, 'lonEdge', lonEdge)
      call mpas_pool_get_array(meshPool, 'latEdgeMid', latEdgeMid)
      call mpas_pool_get_array(meshPool, 'lonEdgeMid', lonEdgeMid)
      call mpas_pool_get_array(meshPool, 'h_s', h_s)
      call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
      call mpas_pool_get_array(meshPool, 'fVertex', fVertex)

      call mpas_pool_get_array(statePool, 'u', u, 1)
      call mpas_pool_get_array(statePool, 'h', h, 1)
      call mpas_pool_get_array(statePool, 'tracers', tracers, 1)
      nTracers = size(tracers, dim=1)

      !HCm = 0.0 !TRiSK
      HCm = 1.0

      xCell = xCell * a
      yCell = yCell * a
      zCell = zCell * a
      xVertex = xVertex * a
      yVertex = yVertex * a
      zVertex = zVertex * a
      xEdge = xEdge * a
      yEdge = yEdge * a
      zEdge = zEdge * a
      xEdgeMid = xEdgeMid * a
      yEdgeMid = yEdgeMid * a
      zEdgeMid = zEdgeMid * a
      xEdgeMidSphere = xEdgeMidSphere * a
      yEdgeMidSphere = yEdgeMidSphere * a
      zEdgeMidSphere = zEdgeMidSphere * a      
      dvEdge = dvEdge * a
      dcEdge = dcEdge * a
      areaCell = areaCell * a**2.0
      areaTriangle = areaTriangle * a**2.0
      kiteAreasOnVertex = kiteAreasOnVertex * a**2.0
      areaPV = areaPV * a**2.0
      areaHEdge = areaHEdge * a**2.0

      !
      ! Initialize wind field
      !
      allocate(psiVertex(nVertices))
      do iVtx = 1, nVertices
         psiVertex(iVtx) = -a * u0 * ( &
                                       sin(latVertex(iVtx)) * cos(alpha) - &
                                       cos(lonVertex(iVtx)) * cos(latVertex(iVtx)) * sin(alpha) &
                                     )
      end do
      do iEdge = 1, nEdges
         u(1,iEdge) = -1.0 * ( &
                                               psiVertex(verticesOnEdge(2,iEdge)) - &
                                               psiVertex(verticesOnEdge(1,iEdge)) &
                                             ) / dvEdge(iEdge)
      end do
      deallocate(psiVertex)

      !
      ! Generate rotated Coriolis field
      !
      if (HCm == 0) then 
         do iEdge = 1, nEdges
            fEdge(iEdge) = 2.0 * omega * &
                                           (-cos(lonEdge(iEdge)) * cos(latEdge(iEdge)) * sin(alpha) + &
                                             sin(latEdge(iEdge)) * cos(alpha) &
                                           )
         end do
      elseif (HCm == 1) then
         do iEdge = 1, nEdges
            fEdge(iEdge) = 2.0 * omega * &
                                           (-cos(lonEdgeMid(iEdge)) * cos(latEdgeMid(iEdge)) * sin(alpha) + &
                                             sin(latEdgeMid(iEdge)) * cos(alpha) &
                                           )
         end do
      end if 
      do iVtx = 1, nVertices
         fVertex(iVtx) = 2.0 * omega * &
                                         (-cos(lonVertex(iVtx)) * cos(latVertex(iVtx)) * sin(alpha) + &
                                           sin(latVertex(iVtx)) * cos(alpha) &
                                         )
      end do

      !
      ! Initialize mountain
      !
      do iCell = 1, nCells
         if (lonCell(iCell) < 0.0) lonCell(iCell) = lonCell(iCell) + 2.0 * pii
         r = sqrt(min(rr**2.0, (lonCell(iCell) - lambda_c)**2.0 + (latCell(iCell) - theta_c)**2.0))
         h_s(iCell) = hs0 * (1.0 - r/rr)
      end do

      !
      ! Initialize tracer fields
      !
      do iCell = 1, nCells
         r = sqrt(min(rr**2.0, (lonCell(iCell) - lambda_c)**2.0 + (latCell(iCell) - theta_c)**2.0))
         tracers(1,1,iCell) = 1.0 - r/rr
      end do
      if (nTracers > 1) then
         do iCell = 1, nCells
            r = sqrt(min(rr**2.0, (lonCell(iCell) - lambda_c)**2.0 + &
                         (latCell(iCell) - theta_c - pii/6.0)**2.0 &
                        ) &
                    )
            tracers(2,1,iCell) = 1.0 - r/rr
         end do
      end if

      !
      ! Initialize height field (actually, fluid thickness field)
      !
      do iCell = 1, nCells
         h(1,iCell) = (gh0 - (a * omega * u0 + 0.5 * u0**2.0) * &
                                         (-cos(lonCell(iCell)) * cos(latCell(iCell)) * sin(alpha) + &
                                           sin(latCell(iCell)) * cos(alpha) &
                                         )**2.0 &
                                      ) / &
                                      gravity
         h(1,iCell) = h(1,iCell) - h_s(iCell)
      end do

   end subroutine sw_test_case_5


   subroutine sw_test_case_6(meshPool, statePool)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup shallow water test case 6: Rossby-Haurwitz Wave
   !
   ! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical 
   !            Approximations to the Shallow Water Equations in Spherical 
   !            Geometry" J. of Comp. Phys., 102, pp. 211--224
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool
      type (mpas_pool_type), intent(inout) :: statePool

      real (kind=RKIND), parameter :: h0 = 8000.0
      real (kind=RKIND), parameter :: w = 7.848e-6
      real (kind=RKIND), parameter :: K = 7.848e-6
      real (kind=RKIND), parameter :: R = 4.0

      integer :: iCell, iEdge, iVtx, nTracers
      integer, pointer :: nVertices, nEdges, nCells
      real (kind=RKIND) :: v
      real (kind=RKIND), allocatable, dimension(:) :: psiVertex

      integer, dimension(:,:), pointer :: verticesOnEdge

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex, xEdge, yEdge, zEdge, dvEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: xEdgeMid, yEdgeMid, zEdgeMid, xEdgeMidSphere, yEdgeMidSphere, zEdgeMidSphere
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex
      real (kind=RKIND), dimension(:), pointer :: areaCell, areaTriangle, h_s
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex

      real (kind=RKIND), dimension(:,:), pointer :: u, h
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers



      !
      ! Scale all distances and areas from a unit sphere to one with radius a
      !

      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
      call mpas_pool_get_array(meshPool, 'xEdgeMid', xEdgeMid)
      call mpas_pool_get_array(meshPool, 'yEdgeMid', yEdgeMid)
      call mpas_pool_get_array(meshPool, 'zEdgeMid', zEdgeMid)
      call mpas_pool_get_array(meshPool, 'xEdgeMidSphere', xEdgeMidSphere)
      call mpas_pool_get_array(meshPool, 'yEdgeMidSphere', yEdgeMidSphere)
      call mpas_pool_get_array(meshPool, 'zEdgeMidSphere', zEdgeMidSphere)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
      call mpas_pool_get_array(meshPool, 'lonVertex', lonVertex)
      call mpas_pool_get_array(meshPool, 'h_s', h_s)

      call mpas_pool_get_array(statePool, 'u', u, 1)
      call mpas_pool_get_array(statePool, 'h', h, 1)
      call mpas_pool_get_array(statePool, 'tracers', tracers, 1)
      nTracers = size(tracers, dim=1)


      xCell = xCell * a
      yCell = yCell * a
      zCell = zCell * a
      xVertex = xVertex * a
      yVertex = yVertex * a
      zVertex = zVertex * a
      xEdgeMid = xEdgeMid * a
      yEdgeMid = yEdgeMid * a
      zEdgeMid = zEdgeMid * a
      xEdgeMidSphere = xEdgeMidSphere * a
      yEdgeMidSphere = yEdgeMidSphere * a
      zEdgeMidSphere = zEdgeMidSphere * a
      xEdge = xEdge * a
      yEdge = yEdge * a
      zEdge = zEdge * a
      dvEdge = dvEdge * a
      dcEdge = dcEdge * a
      areaCell = areaCell * a**2.0
      areaTriangle = areaTriangle * a**2.0
      kiteAreasOnVertex = kiteAreasOnVertex * a**2.0

      !
      ! Initialize wind field
      !
      allocate(psiVertex(nVertices))
      do iVtx = 1, nVertices
         psiVertex(iVtx) = -a * a * w * sin(latVertex(iVtx)) + &
                            a *a * K * (cos(latVertex(iVtx))**R) * &
                            sin(latVertex(iVtx)) * cos(R * lonVertex(iVtx))
      end do
      do iEdge = 1, nEdges
         u(1,iEdge) = -1.0 * ( &
                                               psiVertex(verticesOnEdge(2,iEdge)) - &
                                               psiVertex(verticesOnEdge(1,iEdge)) &
                                             ) / dvEdge(iEdge)
      end do
      deallocate(psiVertex)

      !
      ! Initialize height field (actually, fluid thickness field)
      !
      do iCell = 1, nCells
         h(1,iCell) = (gravity * h0 + a*a*aa(latCell(iCell)) + &
                                                      a*a*bb(latCell(iCell)) * cos(R*lonCell(iCell)) + &
                                                      a*a*cc(latCell(iCell)) * cos(2.0*R*lonCell(iCell)) &
                                      ) / gravity
      end do

   end subroutine sw_test_case_6 

   subroutine sw_test_case_7(meshPool, statePool)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Planar test case with 16 exact hexagons with constant velocity u: to check  
   ! that ke and ke_baricentric are equal when HCm coincides with HCt
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool
      type (mpas_pool_type), intent(inout) :: statePool

      real (kind=RKIND), parameter :: u0 = 2.0 * pii * a / (12.0 * 86400.0)
      real (kind=RKIND), parameter :: gh0 = 29400.0
      real (kind=RKIND), parameter :: alpha = 0.0
      real (kind=RKIND), parameter :: h0 = 1000.0 

      integer :: iCell, cell2, iEdge, jEdge, iVtx, i, j, k, ni, nj 
      integer, pointer :: nVertices, nEdges, nCells, nVertLevels
      real (kind=RKIND) :: v, xVelocity, yVelocity, dotProd, const1, const2, const3  
      real (kind=RKIND), allocatable, dimension(:) :: psiVertex

      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: verticesOnEdge, edgesOnCell, cellsOnEdge, edgeSignOnCell

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex, xEdge, yEdge, zEdge, dvEdge, dcEdge, angleEdge 
      real (kind=RKIND), dimension(:), pointer :: xEdgeMid, yEdgeMid, zEdgeMid, xEdgeMidSphere, yEdgeMidSphere, zEdgeMidSphere
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex
      real (kind=RKIND), dimension(:), pointer :: latEdge, lonEdge
      real (kind=RKIND), dimension(:), pointer :: latEdgeMid, lonEdgeMid
      real (kind=RKIND), dimension(:), pointer :: areaCell, areaTriangle, fEdge, fVertex
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex, ke, ke_baricentric, keInitial

      real (kind=RKIND), dimension(:,:), pointer :: u, h

      !
      ! Scale all distances and areas from a unit sphere to one with radius a
      !

      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
      call mpas_pool_get_array(meshPool, 'xEdgeMid', xEdgeMid)
      call mpas_pool_get_array(meshPool, 'yEdgeMid', yEdgeMid)
      call mpas_pool_get_array(meshPool, 'zEdgeMid', zEdgeMid)
      call mpas_pool_get_array(meshPool, 'xEdgeMidSphere', xEdgeMidSphere)
      call mpas_pool_get_array(meshPool, 'yEdgeMidSphere', yEdgeMidSphere)
      call mpas_pool_get_array(meshPool, 'zEdgeMidSphere', zEdgeMidSphere)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
      call mpas_pool_get_array(meshPool, 'fVertex', fVertex)
      call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
      call mpas_pool_get_array(meshPool, 'lonVertex', lonVertex)
      call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
      call mpas_pool_get_array(meshPool, 'lonEdge', lonEdge)
      call mpas_pool_get_array(meshPool, 'latEdgeMid', latEdgeMid)
      call mpas_pool_get_array(meshPool, 'lonEdgeMid', lonEdgeMid)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)

      call mpas_pool_get_array(statePool, 'u', u, 1)
      call mpas_pool_get_array(statePool, 'h', h, 1)
      call mpas_pool_get_array(statePool, 'ke', ke, 1)
      call mpas_pool_get_array(statePool, 'ke_baricentric', ke_baricentric, 1)
      call mpas_pool_get_array(statePool, 'keInitial', keInitial)     

      !horizontal velocity u
      xVelocity = 1.0  
      yVelocity = 0.0  
      do iEdge = 1, nEdges
         u(1,iEdge) = cos(angleEdge(iEdge))*xVelocity + sin(angleEdge(iEdge))*yVelocity
         !u(1, iEdge) = 0.5 
      end do

      do iCell = 1, nCells
         h(1,iCell) = h0
      end do 

      ke_baricentric(:,:) = 0.0
      ke(:,:) = 0.0   
      do iCell = 1, nCells     
         print*, iCell 
         do i = 1, nEdgesOnCell(iCell)
            const1 = 0.0
            const2 = 0.0
            const3 = 0.0
            iEdge = edgesOnCell(i,iCell)
            ! ---- BEGIN: compute ni 
            if (iCell == cellsOnEdge(1,iEdge)) then
               cell2 = cellsOnEdge(2,iEdge)
            elseif (iCell == cellsOnEdge(2,iEdge)) then
               cell2 = cellsOnEdge(1,iEdge)
            end if
            if (iCell < cell2) then
               ni = 1.0 !this is n_{e',i} in Peixoto
            else if (iCell > cell2) then
               ni = - 1.0
            end if
            ! --- END: compute ni 
            do k = 1, nVertLevels
               print*, 'right ', dcEdge(iEdge)*u(k,iEdge)*0.5 
               !Peixoto 2016 eq. (13)
               ke(k,iCell) = ke(k,iCell) + 0.25 * dcEdge(iEdge) * dvEdge(iEdge) * u(k,iEdge)**2.0 / areaCell(iCell)
               !print*, 'ke(k,iCell)', ke(k,iCell)
            end do 
            do j = 1, nEdgesOnCell(iCell)
               jEdge = edgesOnCell(j,iCell)
               ! ---- BEGIN: compute nj 
               if (iCell == cellsOnEdge(1,jEdge)) then
                  cell2 = cellsOnEdge(2,jEdge)
               elseif (iCell == cellsOnEdge(2,jEdge)) then
                  cell2 = cellsOnEdge(1,jEdge)
               end if
               if (iCell < cell2) then
                  nj = 1.0 !this is n_{e',i} in Peixoto
               else if (iCell > cell2) then
                  nj = - 1.0
               end if 
               ! --- END: compute nj     
               dotProd = ni*(xEdgeMid(iEdge) - xCell(iCell)) * nj*(xEdgeMid(jEdge) - xCell(iCell)) + &
                         ni*(yEdgeMid(iEdge) - yCell(iCell)) * nj*(yEdgeMid(jEdge) - yCell(iCell)) + &
                         ni*(zEdgeMid(iEdge) - zCell(iCell)) * nj*(zEdgeMid(jEdge) - zCell(iCell))
               do k = 1, nVertLevels
                  const1 = const1 + (u(k,jEdge) * dvEdge(jEdge) * nj*(xEdgeMid(jEdge) - xCell(iCell)))/areaCell(iCell)
                  const2 = const2 + (u(k,jEdge) * dvEdge(jEdge) * nj*(yEdgeMid(jEdge) - yCell(iCell)))/areaCell(iCell)
                  const3 = const3 + (u(k,jEdge) * dvEdge(jEdge) * nj*(zEdgeMid(jEdge) - zCell(iCell)))/areaCell(iCell) 
                  !ke_baricentric(k,iCell) = ke_baricentric(k,iCell) + dvEdge(iEdge) * u(k,iEdge) * dvEdge(jEdge) * u(k,jEdge) * dotProd !dvEdge=l, dcEdge=d 
                  ke_baricentric(k,iCell) = ke_baricentric(k,iCell) + 0.5 * (dvEdge(iEdge) * u(k,iEdge) * dvEdge(jEdge) * u(k,jEdge) * dotProd) / areaCell(iCell)**2.0 
                  !print*, 'ke(k,iCell)_baricentric ', ke_baricentric(k,iCell)       
               end do
            end do     
            print*, 'left 1 ', ni*(xEdgeMid(iEdge) - xCell(iCell)) 
            print*, 'left 2 ', const1*ni*(xEdgeMid(iEdge) - xCell(iCell)) + const2*ni*(yEdgeMid(iEdge) - yCell(iCell)) + const3*ni*(zEdgeMid(iEdge) - zCell(iCell))  
            print*, 'diff ', abs(dcEdge(iEdge)*u(1,iEdge)*0.5 - (const1*ni*(xEdgeMid(iEdge) - xCell(iCell)) + const2*ni*(yEdgeMid(iEdge) - yCell(iCell)) + const3*ni*(zEdgeMid(iEdge) - zCell(iCell))))   
         end do
         print*, 'ke ', 'ke_baricentric ', ke(1,iCell), ke_baricentric(1,iCell)
         print*, 'diff ke ', abs(ke_baricentric(1,iCell) - ke(1,iCell))
         print*,'---------------'
      end do

   end subroutine sw_test_case_7

   subroutine ocn_init_routines_setup_sign_and_index_fields(meshPool) 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! This routine initializes edgeSignOnCell, edgeSignOnVertex and
   ! kiteIndexOnCell
   !         
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      type (mpas_pool_type), intent(inout) :: meshPool

      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell, edgesOnVertex, cellsOnVertex, cellsOnEdge, verticesOnCell, verticesOnEdge
      integer, dimension(:,:), pointer :: edgeSignOnCell, edgeSignOnVertex, kiteIndexOnCell

      integer, pointer :: nCells, nEdges, nVertices, vertexDegree
      integer :: iCell, iEdge, iVertex, i, j, k
 
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnVertex', edgesOnVertex)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnVertex', edgeSignOnVertex)
      call mpas_pool_get_array(meshPool, 'kiteIndexOnCell', kiteIndexOnCell)

      edgeSignOnCell = 0.0_RKIND
      edgeSignOnVertex = 0.0_RKIND
      kiteIndexOnCell = 0.0_RKIND

      do iCell = 1, nCells
        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          iVertex = verticesOnCell(i, iCell)

          ! Vector points from cell 1 to cell 2
          if(iCell == cellsOnEdge(1, iEdge)) then
            edgeSignOnCell(i, iCell) = -1
          else
            edgeSignOnCell(i, iCell) =  1
          end if

          do j = 1, vertexDegree
            if(cellsOnVertex(j, iVertex) == iCell) then
              kiteIndexOnCell(i, iCell) = j
            end if
          end do
        end do
      end do

      do iVertex = 1, nVertices
        do i = 1, vertexDegree
          iEdge = edgesOnVertex(i, iVertex)

          ! Vector points from vertex 1 to vertex 2
          if(iVertex == verticesOnEdge(1, iEdge)) then
            edgeSignOnVertex(i, iVertex) = -1
          else
            edgeSignOnVertex(i, iVertex) =  1
          end if
        end do
      end do

   end subroutine ocn_init_routines_setup_sign_and_index_fields

   subroutine compute_weightsOnEdge_fortran(meshPool)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Computing the weights for the Coriolis term in the equation of u 
   !         
   !
   ! Reference: Pedro S. Peixoto, "Accuracy analysis of mimetic finite
   !            volume operators on geodesic grids and a consistent alternative" 
   !            J. of Comp. Phys., 310, pp. 127--160
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool 

      integer :: iEdge, jEdge, i, ii, jj, v2, sharedCell, v2Edge, jEdgeIsNextEdge, nextEdge, nextEdgeOld
      integer, pointer :: nEdges, vertexDegree

      integer, dimension(2) :: iEdgeCells, jEdgeCells, jEdgeVertices, nextEdgeCells
 
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge
      integer, dimension(:,:), pointer :: verticesOnEdge, edgesOnCell, cellsOnEdge, edgesOnEdge, edgesOnVertex

      real (kind=RKIND) :: areaSum, laCell, loCell, laVertex, loVertex, laEdge, loEdge, c, n, zMidpoint, yMidpoint, xMidpoint, deltaLon, Bx, By
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex
      real (kind=RKIND), dimension(:), pointer :: latEdge, lonEdge
      real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge, areaCell, xCell, yCell, zCell

      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdgeFortran
 
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
      call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnVertex', edgesOnVertex)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)

      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
      call mpas_pool_get_array(meshPool, 'lonEdge', lonEdge)
      call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
      call mpas_pool_get_array(meshPool, 'lonVertex', lonVertex)

      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'weightsOnEdgeFortran', weightsOnEdgeFortran)

      weightsOnEdgeFortran(:,:) = 0.0
      do iEdge = 1, nEdges
         iEdgeCells(1) = cellsOnEdge(1,iEdge)
         iEdgeCells(2) = cellsOnEdge(2,iEdge)
         
         v2 = verticesOnEdge(1,iEdge) !this is v2 in the Thuburn paper  
         ! ---- BEGIN: compute c
         do ii = 1, vertexDegree
            v2Edge = edgesOnVertex(ii, v2)
            if (v2Edge == iEdge) then 
               ! Vector points from vertex 1 to vertex 2
               if(v2 == verticesOnEdge(1, iEdge)) then
                  c  = - 1.0 !this is t_{e,v_2} from Thuburn (or edgeSignOnVertex(ii,v2) from MPAS)
               else
                  c  =  1.0
               end if
            end if
         end do       
 
         c = - c !this is c_{e,e'} from Peixoto, which is -1/t_{e,v_2}
         
         ! ----- END: compute c

         do i = 1, nEdgesOnEdge(iEdge)
            areaSum = 0.0
            v2 = verticesOnEdge(1,iEdge)
            jEdge = edgesOnEdge(i,iEdge)
            jEdgeCells(1) = cellsOnEdge(1,jEdge)
            jEdgeCells(2) = cellsOnEdge(2,jEdge)
            do ii = 1,2
               do jj = 1,2
                  if (iEdgeCells(ii)  == jEdgeCells(jj)) then
                     sharedCell = iEdgeCells(ii) !this is to find what is the cell in common (in what cell jEdge is)
                  end if
               end do
            end do         

            ! ---- BEGIN: compute n
            if (sharedCell < jEdgeCells(1) .or. sharedCell < jEdgeCells(2)) then
               n = 1.0 !this is n_{e',i} in Peixoto
            else if (sharedCell > jEdgeCells(1) .or. sharedCell > jEdgeCells(2)) then
               n = - 1.0
            end if
            ! --- END: compute n

            ! BEGIN: compute areaSum 
            ! here we compute the area of the only triangle associated with iEdge
            laCell = latCell(sharedCell)
            loCell = lonCell(sharedCell) 
            laVertex = latVertex(v2)
            loVertex = lonVertex(v2)
               
            deltaLon = lonCell(iEdgeCells(2)) - lonCell(iEdgeCells(1))
            Bx = cos(latCell(iEdgeCells(2))) * cos(deltaLon)
            By = cos(latCell(iEdgeCells(2))) * sin(deltaLon)
            laEdge = atan2(sin(latCell(iEdgeCells(1))) + sin(latCell(iEdgeCells(2))) , sqrt( (cos(latCell(iEdgeCells(1))) + Bx) * (cos(latCell(iEdgeCells(1))) + Bx) + By * By ))
            loEdge = lonCell(iEdgeCells(1)) + atan2(By, cos(latCell(iEdgeCells(1))) + Bx)

            !BEGIN: COMMENT TO DO CELL AREA TEST
            areaSum = areaSum + sphericalTriangleArea(laEdge, loEdge, laVertex, loVertex, laCell, loCell) / areaCell(sharedCell)
            !END: COMMENT TO DO CELL AREA TEST

            !BEGIN: UNCOMMENT TO DO CELL AREA TEST
            !areaSum = areaSum + sphericalTriangleArea(laEdge, loEdge, laVertex, loVertex, laCell, loCell)
            !END: UNCOMMENT TO DO CELL AREA TEST

            jEdgeVertices(1)  = verticesOnEdge(1, jEdge)
            jEdgeVertices(2)  = verticesOnEdge(2, jEdge)
            
            !BEGIN: TO REMOVE
            !print*, 'v2 =', v2, 'jEdgeVertices= ', jEdgeVertices(1), jEdgeVertices(2)
            !END: TO REMOVE
            
            !BEGIN: COMMENT TO DO CELL AREA TEST
            if ((jEdgeVertices(1) == v2) .or. (jEdgeVertices(2) == v2)) then 
               jEdgeIsNextEdge = 1
            else
               jEdgeIsNextEdge = 0
            end if
            !END: COMMENT TO DO CELL AREA TEST

            !BEGIN: UNCOMMENT TO DO CELL AREA TEST
            !jEdgeIsNextEdge = 0
            !END: UNCOMMENT TO DO CELL AREA TEST

            if (jEdgeIsNextEdge == 0) then
               do ii = 1, vertexDegree
                  nextEdge = edgesOnVertex(ii,v2)
                  if ((cellsOnEdge(1,nextEdge) == sharedCell .or. cellsOnEdge(2,nextEdge) == sharedCell) .and. (nextEdge .ne. iEdge )) then
                     exit
                  end if
               end do
            end if

            !BEGIN: TO REMOVE
            !print*, 'BEFORE THE WHILE LOOP'
            !END: TO REMOVE

            nextEdgeOld = iEdge
            do while (jEdgeIsNextEdge == 0)

               !BEGIN: TO REMOVE
               !print*, 'nextEdge =', nextEdge, 'nextEdgeVertices=', verticesOnEdge(1,nextEdge), verticesOnEdge(2,nextEdge) 
               !END: TO REMOVE

               ! this computes the areas of the two triangles associated with nextEdge
               nextEdgeCells(1) = cellsOnEdge(1,nextEdge)
               nextEdgeCells(2) = cellsOnEdge(2,nextEdge)

               deltaLon = lonCell(nextEdgeCells(2)) - lonCell(nextEdgeCells(1))
               Bx = cos(latCell(nextEdgeCells(2))) * cos(deltaLon)
               By = cos(latCell(nextEdgeCells(2))) * sin(deltaLon)
               laEdge = atan2(sin(latCell(nextEdgeCells(1))) + sin(latCell(nextEdgeCells(2))) , sqrt( (cos(latCell(nextEdgeCells(1))) + Bx) * (cos(latCell(nextEdgeCells(1))) + Bx) + By * By ))
               loEdge = lonCell(nextEdgeCells(1)) + atan2(By, cos(latCell(nextEdgeCells(1))) + Bx)

               do ii=1,2 
                  laVertex = latVertex(verticesOnEdge(ii, nextEdge))
                  loVertex = lonVertex(verticesOnEdge(ii, nextEdge))

                  !BEGIN: COMMENT TO DO CELL AREA TEST
                  areaSum = areaSum + sphericalTriangleArea(laEdge, loEdge, laVertex, loVertex, laCell, loCell) / areaCell(sharedCell)
                  !END: COMMENT TO DO CELL AREA TEST

                  !BEGIN: UNCOMMENT TO DO CELL AREA TEST
                  !areaSum = areaSum + sphericalTriangleArea(laEdge, loEdge, laVertex, loVertex, laCell, loCell)
                  !END: UNCOMMENT TO DO CELL AREA TEST
             
               end do

               ! new we have to find the nextEdge
               do ii=1,2
                  if (verticesOnEdge(ii,nextEdge) .ne. v2) then
                     v2 = verticesOnEdge(ii,nextEdge) !this is necessary to find the nextEdge
                     !BEGIN: TO REMOVE
                     !print*, 'v2=', v2
                     !END: TO REMOVE
                     exit
                  end if
               end do

               !BEGIN: TO REMOVE
               !print*, 'iEdgeVertices =', v2, verticesOnEdge(2,iEdge), 'jEdgeVertices= ', jEdgeVertices(1), jEdgeVertices(2)
               !END: TO REMOVE

               !BEGIN: COMMENT TO DO CELL AREA TEST
               if (jEdgeVertices(1) == v2 .or. jEdgeVertices(2) == v2) then
                  jEdgeIsNextEdge = 1
               end if
               !END: COMMENT TO DO CELL AREA TEST

               !BEGIN: UNCOMMENT TO DO CELL AREA TEST
               !if (v2 == verticesOnEdge(2,iEdge)) then
               !   jEdgeIsNextEdge = 1
               !end if
               !END: UNCOMMENT TO DO CELL AREA TEST

               if (jEdgeIsNextEdge == 0) then
                  nextEdgeOld = nextEdge
                  do ii = 1, vertexDegree
                     nextEdge = edgesOnVertex(ii,v2)
                     if ((cellsOnEdge(1,nextEdge) == sharedCell .or. cellsOnEdge(2,nextEdge) == sharedCell) .and. (nextEdge .ne. nextEdgeOld )) then
                        exit
                     end if
                  end do
               end if

            end do
 
            !BEGIN: TO REMOVE
            !print*, 'AFTER THE WHILE LOOP'
            !END: TO REMOVE
            
            ! here we compute the area of the only triangle associated with jEdge
         
            laVertex = latVertex(v2)
            loVertex = lonVertex(v2)

            !BEGIN: COMMENT TO DO CELL AREA TEST
            deltaLon = lonCell(jEdgeCells(2)) - lonCell(jEdgeCells(1))
            Bx = cos(latCell(jEdgeCells(2))) * cos(deltaLon)
            By = cos(latCell(jEdgeCells(2))) * sin(deltaLon)
            laEdge = atan2(sin(latCell(jEdgeCells(1))) + sin(latCell(jEdgeCells(2))) , sqrt( (cos(latCell(jEdgeCells(1))) + Bx) * (cos(latCell(jEdgeCells(1))) + Bx) + By * By ))
            loEdge = lonCell(jEdgeCells(1)) + atan2(By, cos(latCell(jEdgeCells(1))) + Bx)            
            !END: COMMENT TO DO CELL AREA TEST

            !BEGIN: UNCOMMENT TO DO CELL AREA TEST
            !deltaLon = lonCell(iEdgeCells(2)) - lonCell(iEdgeCells(1))
            !Bx = cos(latCell(iEdgeCells(2))) * cos(deltaLon)
            !By = cos(latCell(iEdgeCells(2))) * sin(deltaLon)
            !laEdge = atan2(sin(latCell(iEdgeCells(1))) + sin(latCell(iEdgeCells(2))) , sqrt( (cos(latCell(iEdgeCells(1))) + Bx) * (cos(latCell(iEdgeCells(1))) + Bx) + By * By ))
            !loEdge = lonCell(iEdgeCells(1)) + atan2(By, cos(latCell(iEdgeCells(1))) + Bx)            
            !END: UNCOMMENT TO DO CELL AREA TEST

            !BEGIN: COMMENT TO DO CELL AREA TEST
            areaSum = areaSum + sphericalTriangleArea(laEdge, loEdge, laVertex, loVertex, laCell, loCell) / areaCell(sharedCell)
            !END: COMMENT TO DO CELL AREA TEST

            !BEGIN: UNCOMMENT TO DO CELL AREA TEST
            !areaSum = areaSum + sphericalTriangleArea(laEdge, loEdge, laVertex, loVertex, laCell, loCell)
            !print*, 'areaSum =', areaSum, 'areaCell =', areaCell(sharedCell), 'difference =', abs(areaSum - areaCell(sharedCell)), 'edge =', iEdge, 'cell =', sharedCell
            !END: UNCOMMENT TO DO CELL AREA TEST

            ! --- END: compute areaSum

            weightsOnEdgeFortran(i, iEdge) = c  * dvEdge(jEdge)/dcEdge(iEdge) * (0.5 - areaSum) * n 

         end do
      end do


   end subroutine compute_weightsOnEdge_fortran

   subroutine compute_weightsOnEdge_baricentric(meshPool)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Computing the weights for the Coriolis term in the equation of u using eq.
   ! (37) in Peixoto. Assumptions made: 
   ! 1. t_e = (x_e - x_v)/||x_e - x_v||
   ! 2. we assume that n_e goes from the center of cellsOnEdge(1, iEdge) to
   ! cellsOnEdge(2, iEdge)  
   !         
   !
   ! Reference: Pedro S. Peixoto, "Accuracy analysis of mimetic finite
   !            volume operators on geodesic grids and a consistent alternative" 
   !            J. of Comp. Phys., 310, pp. 127--160
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool

      integer :: iEdge, jEdge, i, j, ii, jj, v2, sharedCell, v2Edge, jEdgeIsNextEdge, nextEdge, nextEdgeOld, iVertex, n, indexEdge  
      integer, pointer :: nEdges, vertexDegree

      integer, dimension(2) :: iEdgeCells, jEdgeCells, jEdgeVertices, nextEdgeCells

      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge, nEdgesArray
      integer, dimension(:,:), pointer :: verticesOnEdge, edgesOnCell, cellsOnEdge, edgesOnEdge, edgesOnVertex, edgeSignOnCell 

      real (kind=RKIND) :: zMidpoint, yMidpoint, xMidpoint, xTangent, yTangent, zTangent, normTangent, zCross, xCross, yCross, dotProd, ne_x, ne_y, ne_z  
      real (kind=RKIND), dimension(:), pointer :: xEdgeMidSphere, yEdgeMidSphere, zEdgeMidSphere, xEdge, yEdge, zEdge
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
      real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge, areaCell, xCell, yCell, zCell

      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge_baricentric

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray) 
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
      call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnVertex', edgesOnVertex)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)

      call mpas_pool_get_array(meshPool, 'xEdgeMidSphere', xEdgeMidSphere)
      call mpas_pool_get_array(meshPool, 'yEdgeMidSphere', yEdgeMidSphere)
      call mpas_pool_get_array(meshPool, 'zEdgeMidSphere', zEdgeMidSphere)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(meshPool, 'zEdge', zEdge)      
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)

      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'weightsOnEdge_baricentric', weightsOnEdge_baricentric)

      weightsOnEdge_baricentric(:,:) = 0.0
      !print*, 'nEdgesArray ', nEdgesArray
      !print*, 'nEdges ', nEdges  
      !do iEdge = 1, nEdges
      do iEdge = 1, nEdgesArray(3) 
         iEdgeCells(1) = cellsOnEdge(1,iEdge)
         iEdgeCells(2) = cellsOnEdge(2,iEdge)
         !--- begin t_e Darren ---
         !do j = 1, nEdgesOnCell(iEdgeCells(1))
         !   if (iEdge == edgesOnCell(j,iEdgeCells(1))) then
         !      indexEdge = j
         !   end if
         !end do
         ! ---- BEGIN: compute n
         if (iEdgeCells(1) > iEdgeCells(2)) then 
            n = - 1.0 
         elseif (iEdgeCells(2) > iEdgeCells(1)) then 
            n = 1.0 
         end if 
         ! ---- END: compute n
         ne_x = (xCell(iEdgeCells(1)) - xCell(iEdgeCells(2)))
         ne_y = (yCell(iEdgeCells(1)) - yCell(iEdgeCells(2)))
         ne_z = (zCell(iEdgeCells(1)) - zCell(iEdgeCells(2)))
         xTangent = yEdgeMidSphere(iEdge) * ne_z - zEdgeMidSphere(iEdge) * ne_y  
         yTangent = zEdgeMidSphere(iEdge) * ne_x - xEdgeMidSphere(iEdge) * ne_z 
         zTangent = xEdgeMidSphere(iEdge) * ne_y - yEdgeMidSphere(iEdge) * ne_x
         normTangent = sqrt(xTangent**2.0 + yTangent**2.0 + zTangent**2.0)     
         xTangent = xTangent/normTangent
         yTangent = yTangent/normTangent 
         zTangent = zTangent/normTangent   
         !print*, xTangent, yTangent, zTangent
         !--- end t_e Darren ---
         !--- begin t_e --- 
         !iVertex = verticesOnEdge(1, iEdge)
         !xTangent = xVertex(iVertex) - xEdgeMidSphere(iEdge)
         !yTangent = yVertex(iVertex) - yEdgeMidSphere(iEdge)
         !zTangent = zVertex(iVertex) - zEdgeMidSphere(iEdge) 
         !normTangent = sqrt(xTangent**2.0 + yTangent**2.0 + zTangent**2.0)
         !xTangent = xTangent/normTangent 
         !yTangent = yTangent/normTangent 
         !zTangent = zTangent/normTangent 
         !zCross = xTangent * (yEdgeMidSphere(iEdge) - yCell(iEdgeCells(1))) - yTangent * (xEdgeMidSphere(iEdge) - xCell(iEdgeCells(1)))
         !xCross = yTangent * (zEdgeMidSphere(iEdge) - zCell(iEdgeCells(1))) - zTangent * (yEdgeMidSphere(iEdge) - yCell(iEdgeCells(1)))
         !yCross = zTangent * (xEdgeMidSphere(iEdge) - xCell(iEdgeCells(1))) - xTangent * (zEdgeMidSphere(iEdge) - zCell(iEdgeCells(1)))
         !print*, 'xCross= ', xCross, 'yCross= ', yCross
         !dotProd = xCross * xTangent + &
                   !yCross * yTangent + &
                   ! zCross * zTangent 
         !print*, 'dotProd', dotProd
         !dotProd = (xEdgeMidSphere(iEdge) - xCell(iEdgeCells(1))) * xTangent + &
                   !(yEdgeMidSphere(iEdge) - yCell(iEdgeCells(1))) * yTangent + &
                   !(zEdgeMidSphere(iEdge) - zCell(iEdgeCells(1))) * zTangent
         !print*, 'dotProd', dotProd
         !if (zCross < 0.0) then
         !   xTangent = -1.0 * xTangent
         !   yTangent = -1.0 * yTangent
         !   zTangent = -1.0 * zTangent
         !end if  
         !--- end t_e ---  
         do i = 1, nEdgesOnEdge(iEdge) 
            jEdge = edgesOnEdge(i,iEdge)         
            jEdgeCells(1) = cellsOnEdge(1,jEdge)
            jEdgeCells(2) = cellsOnEdge(2,jEdge)
            do ii = 1,2
               do jj = 1,2
                  if (iEdgeCells(ii)  == jEdgeCells(jj)) then
                     sharedCell = iEdgeCells(ii) !this is to find what is the cell in common (in what cell jEdge is)
                  end if
               end do
            end do
            ! ---- BEGIN: compute n
            if (sharedCell < jEdgeCells(1) .or. sharedCell < jEdgeCells(2)) then
               n = 1.0 !this is n_{e',i} in Peixoto
            else if (sharedCell > jEdgeCells(1) .or. sharedCell > jEdgeCells(2)) then
               n = - 1.0
            end if
            ! --- END: compute n 
            do j = 1, nEdgesOnCell(sharedCell)
               if (jEdge == edgesOnCell(j,sharedCell)) then 
                  indexEdge = j
               end if 
            end do              
            dotProd = edgeSignOnCell(indexEdge, sharedCell)*(xEdgeMidSphere(jEdge) - xCell(sharedCell)) * xTangent + & 
                      edgeSignOnCell(indexEdge, sharedCell)*(yEdgeMidSphere(jEdge) - yCell(sharedCell)) * yTangent + & 
                      edgeSignOnCell(indexEdge, sharedCell)*(zEdgeMidSphere(jEdge) - zCell(sharedCell)) * zTangent 
            !print*, 'dotProd', dotProd
            !print*, 'areaCell(sharedCell) ', areaCell(sharedCell)  
            !print*, 'dvEdge(jEdge)/areaCell(sharedCell)', dvEdge(jEdge)/areaCell(sharedCell) 
            weightsOnEdge_baricentric(i, iEdge) = 0.5 * (dvEdge(jEdge)/areaCell(sharedCell)) * dotProd 
            !print*, weightsOnEdge_baricentric(i, iEdge)  
         end do 
      end do 
   
   end subroutine compute_weightsOnEdge_baricentric

   subroutine xyzEdgeMid(meshPool)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Computation of xEdgeMid, yEdgeMid and zEdgeMid, i.e. computation of the 
   ! x-coordinate, y-coordinate and z-coordinate of the midpoint of a Voronoi edge 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool
      integer, pointer :: nEdges 
      integer :: iEdge, vertex1, vertex2
      real (kind=RKIND) :: xVertex1, yVertex1, zVertex1, xVertex2, yVertex2, zVertex2, length, radius
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex    
      real (kind=RKIND), dimension(:), pointer :: xEdgeMid, yEdgeMid, zEdgeMid, xEdgeMidSphere, yEdgeMidSphere, zEdgeMidSphere 
      integer, dimension(:,:), pointer :: verticesOnEdge

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'xEdgeMid', xEdgeMid)
      call mpas_pool_get_array(meshPool, 'yEdgeMid', yEdgeMid)
      call mpas_pool_get_array(meshPool, 'zEdgeMid', zEdgeMid)
      call mpas_pool_get_array(meshPool, 'xEdgeMidSphere', xEdgeMidSphere)
      call mpas_pool_get_array(meshPool, 'yEdgeMidSphere', yEdgeMidSphere)
      call mpas_pool_get_array(meshPool, 'zEdgeMidSphere', zEdgeMidSphere)

      xEdgeMid(:) = 0.0
      yEdgeMid(:) = 0.0
      zEdgeMid(:) = 0.0
      xEdgeMidSphere(:) = 0.0
      yEdgeMidSphere(:) = 0.0
      zEdgeMidSphere(:) = 0.0
      radius = 6371220.0

      do iEdge = 1, nEdges

         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)  
         xVertex1 = xVertex(vertex1)
         yVertex1 = yVertex(vertex1)
         zVertex1 = zVertex(vertex1)
         xVertex2 = xVertex(vertex2)
         yVertex2 = yVertex(vertex2)
         zVertex2 = zVertex(vertex2)       
         xEdgeMid(iEdge) = 0.5*(xVertex1 + xVertex2) 
         yEdgeMid(iEdge) = 0.5*(yVertex1 + yVertex2)
         zEdgeMid(iEdge) = 0.5*(zVertex1 + zVertex2)
         length = sqrt(xEdgeMid(iEdge)**2.0 + yEdgeMid(iEdge)**2.0 + zEdgeMid(iEdge)**2.0)
         xEdgeMidSphere(iEdge) = (xEdgeMid(iEdge)/length)*radius
         yEdgeMidSphere(iEdge) = (yEdgeMid(iEdge)/length)*radius
         zEdgeMidSphere(iEdge) = (zEdgeMid(iEdge)/length)*radius
         
      end do 

   end subroutine xyzEdgeMid


   subroutine computeLatLonEdgeMid(meshPool)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! computeLatLonEdgeMid computes the latitudine and longitude of the midpoint of 
   ! a Voronoi edge
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool
      integer, pointer :: nEdges 
      integer :: iEdge, vertex1, vertex2 
      real (kind=RKIND) :: laVertex1, loVertex1, laVertex2, loVertex2, deltaLon, Bx, By
      real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex 
      real (kind=RKIND), dimension(:), pointer :: latEdgeMid, lonEdgeMid 
      integer, dimension(:,:), pointer :: verticesOnEdge 

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges) 
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge) 
      call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
      call mpas_pool_get_array(meshPool, 'lonVertex', lonVertex)
      call mpas_pool_get_array(meshPool, 'latEdgeMid', latEdgeMid)
      call mpas_pool_get_array(meshPool, 'lonEdgeMid', lonEdgeMid)

      latEdgeMid(:) = 0.0
      lonEdgeMid(:) = 0.0

      do iEdge = 1, nEdges

         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)      
         laVertex1 = latVertex(vertex1)
         loVertex1 = lonVertex(vertex1)
         laVertex2 = latVertex(vertex2)
         loVertex2 = lonVertex(vertex2)
         deltaLon = loVertex2 - loVertex1
         Bx = cos(laVertex2) * cos(deltaLon)
         By = cos(laVertex2) * sin(deltaLon)
         latEdgeMid(iEdge) = atan2(sin(laVertex1) + sin(laVertex2) , sqrt((cos(laVertex1) + Bx) * (cos(laVertex1) + Bx) + By * By ))
         lonEdgeMid(iEdge) = loVertex1 + atan2(By, cos(laVertex1) + Bx) 
   
      end do

   end subroutine computeLatLonEdgeMid 

   subroutine weightsPV(meshPool)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! weightsPV computes the weights needed in the HCm discretization
   ! for the computation of h_vertex_baricentric 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool
      integer :: iVertex, cell1, cell2, cell3 
      real (kind=RKIND) :: laVertex1, loVertex1, laCell1, loCell1, laCell2,  loCell2, laCell3, loCell3
      real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      integer, dimension(:), pointer :: nVerticesArray 
      integer, dimension(:,:), pointer :: cellsOnVertex
      real (kind=RKIND), dimension(:,:), pointer :: areaPV 
 
      call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
      call mpas_pool_get_array(meshPool, 'lonVertex', lonVertex)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_dimension(meshPool, 'nVerticesArray', nVerticesArray)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'areaPV', areaPV)

      areaPV(:,:) = 0.0

      do iVertex = 1, nVerticesArray(3)
         cell1 = cellsOnVertex(1,iVertex)
         cell2 = cellsOnVertex(2,iVertex)
         cell3 = cellsOnVertex(3,iVertex)
         laCell1 = latCell(cell1)
         loCell1 = lonCell(cell1)
         laCell2 = latCell(cell2)
         loCell2 = lonCell(cell2)
         laCell3 = latCell(cell3)
         loCell3 = lonCell(cell3)
         laVertex1 = latVertex(iVertex)
         loVertex1 = lonVertex(iVertex)
         areaPV(1, iVertex) = sphericalTriangleArea(laCell2, loCell2, laCell3, loCell3, laVertex1, loVertex1)
         areaPV(2, iVertex) = sphericalTriangleArea(laCell1, loCell1, laCell3, loCell3, laVertex1, loVertex1)
         areaPV(3, iVertex) = sphericalTriangleArea(laCell1, loCell1, laCell2, loCell2, laVertex1, loVertex1)
         !areaPV(iVertex, cell1) = sphericalTriangleArea(laCell2, loCell2, laCell3, loCell3, laVertex1, loVertex1)
         !areaPV(iVertex, cell2) = sphericalTriangleArea(laCell1, loCell1, laCell3, loCell3, laVertex1, loVertex1)
         !areaPV(iVertex, cell3) = sphericalTriangleArea(laCell1, loCell1, laCell2, loCell2, laVertex1, loVertex1)
         !area3 = areaTriangle(iVertex) - area1 - area2
      end do 

   end subroutine weightsPV

   subroutine weightsHEDGE(meshPool)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! weightsHEDGE computes the weights needed in the HCm discretization
   ! for the computation of h_edge_baricentric 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool
      integer :: vertex1, vertex2, cell1, cell2, cell3, cellNew, iEdge, triangleCenter, i, k
      real (kind=RKIND) :: laVertex1, loVertex1, laVertex2, loVertex2, laEdge, loEdge, laEdgeMidVoronoi, loEdgeMidVoronoi, & 
                           laCell1, loCell1, laCell2, loCell2, laCell3, loCell3, a, b, radius1
      real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex, latCell, lonCell
      real (kind=RKIND), dimension(:), pointer :: latEdge, lonEdge, latEdgeMid, lonEdgeMid  
      integer, pointer :: nVertLevels, vertexDegree
      integer, dimension(:), pointer :: nEdgesArray
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, verticesOnEdge 
      real (kind=RKIND), dimension(:,:), pointer :: areaHEdge
      real (kind=RKIND), dimension(:), pointer :: areaTriangle
      integer, dimension(:), pointer :: cell3HEdge, triangleHEdge  

      call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
      call mpas_pool_get_array(meshPool, 'lonVertex', lonVertex)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
      call mpas_pool_get_array(meshPool, 'lonEdge', lonEdge)
      call mpas_pool_get_array(meshPool, 'latEdgeMid', latEdgeMid)
      call mpas_pool_get_array(meshPool, 'lonEdgeMid', lonEdgeMid)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'areaHEdge', areaHEdge)
      call mpas_pool_get_array(meshPool, 'cell3HEdge', cell3HEdge)
      call mpas_pool_get_array(meshPool, 'triangleHEdge', triangleHEdge)
      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree) 

      areaHEdge(:,:) = 0.0
      cell3HEdge(:) = 0.0
      triangleHEdge(:) = 0.0

      radius1 = 6371220.0

      do iEdge = 1, nEdgesArray(3)
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge) 
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)
         laVertex1 = latVertex(vertex1)
         loVertex1 = lonVertex(vertex1)
         laVertex2 = latVertex(vertex2)
         loVertex2 = lonVertex(vertex2)
         laEdge = latEdge(iEdge)
         loEdge = lonEdge(iEdge)
         !deltaLon = lonCell(cell2) - lonCell(cell1)
         !Bx = cos(latCell(cell2)) * cos(deltaLon)
         !By = cos(latCell(cell2)) * sin(deltaLon)
         !laEdgeMidTriangle = atan2(sin(latCell(cell1)) +
         !sin(latCell(cell2)) , sqrt( (cos(latCell(cell1)) + Bx) *
         !(cos(latCell(cell1)) + Bx) + By * By ))
         !loEdgeMidTriangle = lonCell(cell1) +
         !atan2(By,cos(latCell(cell1)) + Bx)        
         !print*, 'print mid point triangle and intersection MPAS',
         !laEdgeMidTriangle, loEdgeMidTriangle, laEdge, loEdge 
         !print*, 'difference mid point triangle and intersection MPAS',
         !abs(laEdgeMidTriangle-laEdge), abs(loEdgeMidTriangle-loEdge)
         a = sphere_distance(laVertex1, loVertex1, laEdge, loEdge, radius1)
         b = sphere_distance(laVertex2, loVertex2, laEdge, loEdge, radius1) 
         !print*, 'a', 'b', 'abs(a-b)', a, b, abs(a-b)
         !if (abs(a-b)<1.E-9) then
         if (abs(a-b)<1.E-2) then
            do k = 1, nVertLevels
               !h_edge_baricentric(k,iEdge) = h_edge(k,iEdge)
               areaHEdge(1,iEdge) = 0.5
               areaHEdge(2,iEdge) = 0.5
               areaHEdge(3,iEdge) = 0.0 
            end do
            !print*, 'mid triangle = mid voronoi' 
         else
            if (a<b) then
               triangleCenter = vertex2
            else
               triangleCenter = vertex1
            end if
            do i = 1, vertexDegree
               cellNew = cellsOnVertex(i, triangleCenter)
               if ((cell1.NE.cellNew) .AND. (cell2.NE.cellNew)) then
                  cell3 = cellNew
                  cell3HEdge(iEdge) = cell3 
               end if
            end do
            !print*, 'cell1, cell2, cell3', cell1, cell2, cell3 
            triangleHEdge(iEdge) = triangleCenter
            laEdgeMidVoronoi = latEdgeMid(iEdge)
            loEdgeMidVoronoi = lonEdgeMid(iEdge) 
            !print*, 'print mid point voronoi and triangle',
            !laEdgeMidVoronoi, laEdge, loEdgeMidVoronoi, loEdge
            !print*, 'difference mid point voronoi e triangle',
            !abs(laEdgeMidVoronoi-laEdge), abs(loEdgeMidVoronoi-loEdge) 
            laCell1 = latCell(cell1)
            loCell1 = lonCell(cell1)
            laCell2 = latCell(cell2)
            loCell2 = lonCell(cell2)
            laCell3 = latCell(cell3)
            loCell3 = lonCell(cell3)
            !print*, areaTriangle(triangleCenter)
            areaHEdge(1,iEdge) = sphericalTriangleArea(laCell2, loCell2, laCell3, loCell3, laEdgeMidVoronoi, loEdgeMidVoronoi)
            !print*, areaHEdge(1,iEdge)
            areaHEdge(2,iEdge) = sphericalTriangleArea(laCell1, loCell1, laCell3, loCell3, laEdgeMidVoronoi, loEdgeMidVoronoi)
            !print*, areaHEdge(2,iEdge)
            areaHEdge(3,iEdge) = sphericalTriangleArea(laCell1, loCell1, laCell2, loCell2, laEdgeMidVoronoi, loEdgeMidVoronoi)  
            !print*, areaHEdge(3,iEdge) 
            !print*, '----------'
         end if 
      end do 

   end subroutine weightsHEDGE

   real function sphericalTriangleArea(laEdge, loEdge, laVertex, loVertex, laCell, loCell)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! sphericalTriangleArea uses the spherical analog of Heron's formula to compute
   ! the area of a triangle on the surface of a sphere
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      implicit none
 
      real (kind=RKIND), intent(in) :: laEdge, loEdge, laVertex, loVertex, laCell, loCell
      real (kind=RKIND) :: tanqe, s, a, b, c, one, radius
      one = 1.0 
      radius = 6371220.0  
       
      a = sphere_distance(laCell, loCell, laVertex, loVertex, one)
      b = sphere_distance(laCell, loCell, laEdge, loEdge, one)
      c = sphere_distance(laEdge, loEdge, laVertex, loVertex, one)
      s = 0.5*(a+b+c) 
      !print*, 'printing a, b, c e s', a, b, c, s 
      !print*, 'printing tan', tan(0.5*s), tan(0.5*(s-a)), tan(0.5*(s-b)), tan(0.5*(s-c))

      tanqe = sqrt(tan(0.5*s)*tan(0.5*(s-a))*tan(0.5*(s-b))*tan(0.5*(s-c))) 
      sphericalTriangleArea = 4.*atan(tanqe)*radius**2
      
      !print*, 'printing a, b, c, s e tanqe', a, b, c, s, tanqe

   end function sphericalTriangleArea 

   real function sphere_distance(lat1, lon1, lat2, lon2, radius)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) on a
   !   sphere with given radius.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      implicit none

      real (kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2, radius

      real (kind=RKIND) :: arg1

      arg1 = sqrt( sin(0.5*(lat2-lat1))**2 +  &
                   cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
      sphere_distance = 2.*radius*asin(arg1)

   end function sphere_distance


   real function aa(theta)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! A, used in height field computation for Rossby-Haurwitz wave
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), parameter :: w = 7.848e-6
      real (kind=RKIND), parameter :: K = 7.848e-6
      real (kind=RKIND), parameter :: R = 4.0

      real (kind=RKIND), intent(in) :: theta

      aa = 0.5 * w * (2.0 * omega + w) * cos(theta)**2.0 + &
          0.25 * K**2.0 * cos(theta)**(2.0*R) * ((R+1.0)*cos(theta)**2.0 + 2.0*R**2.0 - R - 2.0 - 2.0*R**2.0 * cos(theta)**(-2.0))

   end function aa

   
   real function bb(theta)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! B, used in height field computation for Rossby-Haurwitz wave
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), parameter :: w = 7.848e-6
      real (kind=RKIND), parameter :: K = 7.848e-6
      real (kind=RKIND), parameter :: R = 4.0

      real (kind=RKIND), intent(in) :: theta

      bb = (2.0*(omega + w)*K / ((R+1.0)*(R+2.0))) * cos(theta)**R * ((R**2.0 + 2.0*R + 2.0) - ((R+1.0)*cos(theta))**2.0)

   end function bb


   real function cc(theta)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C, used in height field computation for Rossby-Haurwitz wave
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), parameter :: w = 7.848e-6
      real (kind=RKIND), parameter :: K = 7.848e-6
      real (kind=RKIND), parameter :: R = 4.0

      real (kind=RKIND), intent(in) :: theta

      cc = 0.25 * K**2.0 * cos(theta)**(2.0*R) * ((R+1.0)*cos(theta)**2.0 - R - 2.0)

   end function cc

end module sw_test_cases
