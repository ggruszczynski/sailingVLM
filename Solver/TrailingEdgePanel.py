from Solver.Panel import Panel
from Solver.vortices import v_induced_by_horseshoe_vortex_basic, v_induced_by_horseshoe_vortex_improved


class TrailingEdgePanel(Panel):
    def __init__(self, p1, p2, p3, p4, gamma_orientation=1):
        super().__init__(p1, p2, p3, p4, gamma_orientation=gamma_orientation)

    def get_induced_velocity(self, ctr_p, V_app_infw):
        [A, B, C, D] = self.get_vortex_ring_position()
        v = v_induced_by_horseshoe_vortex_basic(ctr_p, B, C, V_app_infw, self.gamma_orientation)
        # v = v_induced_by_horseshoe_vortex_improved(ctr_p, A, B, C, D, V_app_infw, self.gamma_orientation)
        return v

    def get_vortex_ring_position(self):
        """
        For a given panel defined by points P1, P2, P3 and P4
        returns the position of the horseshoe vortex defined
        by points A, B and its control point P.

                  ^
                 y|                Points defining the panel
                  |                are named clockwise.
         P3--C----|---P4---D=P3next
          |  |    |    |   |
          |  |    |    |   |
          |  |    +----|---------->
          |  |         |   |      x
          |  |         |   |
         P2--B---------P1--A=P2next

        Parameters
        ----------
        P1, P2, P3, P4 : array_like
                         Points that define the panel

        Returns
        -------
        results : dict
            A, B, C, D - points that define the vortex ring
        """
        p2_p1 = self.p1 - self.p2
        p3_p4 = self.p4 - self.p3

        A = self.p1 + p2_p1 / 4.
        B = self.p2 + p2_p1 / 4.
        C = self.p3 + p3_p4 / 4.
        D = self.p4 + p3_p4 / 4.
        return [A, B, C, D]
