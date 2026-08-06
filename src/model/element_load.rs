//! Element load definitions for 2D finite element models.

use crate::FemError;

/// Coordinate system used to interpret an element load vector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LoadCoordinateSystem2D {
    /// Load components are defined in the model's global x/y coordinate system.
    Global,

    /// Load components are defined in the load-specific local x/y coordinate system.
    Local,
}

/// Uniform line load applied along a 2D beam element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BeamUniformLineLoad2D {
    element_id: usize,
    coordinate_system: LoadCoordinateSystem2D,
    x_component: f64,
    y_component: f64,
}

impl BeamUniformLineLoad2D {
    /// Creates a uniform beam line load.
    ///
    /// The components are force per unit beam length. In the local coordinate
    /// system, `x_component` acts along the beam axis and `y_component` acts
    /// in the local transverse direction.
    pub fn new(
        element_id: usize, coordinate_system: LoadCoordinateSystem2D, x_component: f64, y_component: f64,
    ) -> Result<Self, FemError> {
        if !x_component.is_finite() {
            return Err(FemError::InvalidElementLoadValue {
                element_id,
                load_type: Self::LOAD_TYPE,
                component: "x component",
                value: x_component,
            });
        }

        if !y_component.is_finite() {
            return Err(FemError::InvalidElementLoadValue {
                element_id,
                load_type: Self::LOAD_TYPE,
                component: "y component",
                value: y_component,
            });
        }

        Ok(Self { element_id, coordinate_system, x_component, y_component })
    }

    pub(crate) const LOAD_TYPE: &'static str = "beam_uniform_line";

    /// Returns the loaded element ID.
    #[must_use]
    pub fn element_id(&self) -> usize {
        self.element_id
    }

    /// Returns the coordinate system used by the load components.
    #[must_use]
    pub fn coordinate_system(&self) -> LoadCoordinateSystem2D {
        self.coordinate_system
    }

    /// Returns the x component of the line load.
    #[must_use]
    pub fn x_component(&self) -> f64 {
        self.x_component
    }

    /// Returns the y component of the line load.
    #[must_use]
    pub fn y_component(&self) -> f64 {
        self.y_component
    }

    pub(crate) fn local_components(&self, cosine: f64, sine: f64) -> (f64, f64) {
        match self.coordinate_system {
            LoadCoordinateSystem2D::Local => (self.x_component, self.y_component),
            LoadCoordinateSystem2D::Global => (
                cosine * self.x_component + sine * self.y_component,
                -sine * self.x_component + cosine * self.y_component,
            ),
        }
    }

    pub(crate) fn local_equivalent_nodal_load(&self, length: f64, cosine: f64, sine: f64) -> [f64; 6] {
        let (x_component, y_component) = self.local_components(cosine, sine);

        [
            x_component * length / 2.0,
            y_component * length / 2.0,
            y_component * length.powi(2) / 12.0,
            x_component * length / 2.0,
            y_component * length / 2.0,
            -y_component * length.powi(2) / 12.0,
        ]
    }
}

/// Uniform traction applied along one edge of a plane-stress T3 element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgeTraction2D {
    element_id: usize,
    edge_node_ids: [usize; 2],
    coordinate_system: LoadCoordinateSystem2D,
    x_component: f64,
    y_component: f64,
}

impl EdgeTraction2D {
    /// Creates a uniform edge traction.
    ///
    /// The edge is identified by the two node IDs at its ends. In the local
    /// coordinate system, `x_component` acts from the first edge node toward
    /// the second edge node, and `y_component` acts normal to that edge.
    pub fn new(
        element_id: usize, edge_node_ids: [usize; 2], coordinate_system: LoadCoordinateSystem2D, x_component: f64,
        y_component: f64,
    ) -> Result<Self, FemError> {
        if edge_node_ids[0] == edge_node_ids[1] {
            return Err(FemError::InvalidElementLoadEdge {
                element_id,
                load_type: Self::LOAD_TYPE,
                node_ids: edge_node_ids.to_vec(),
                expected: "two distinct edge node IDs",
            });
        }

        if !x_component.is_finite() {
            return Err(FemError::InvalidElementLoadValue {
                element_id,
                load_type: Self::LOAD_TYPE,
                component: "x component",
                value: x_component,
            });
        }

        if !y_component.is_finite() {
            return Err(FemError::InvalidElementLoadValue {
                element_id,
                load_type: Self::LOAD_TYPE,
                component: "y component",
                value: y_component,
            });
        }

        Ok(Self { element_id, edge_node_ids, coordinate_system, x_component, y_component })
    }

    pub(crate) const LOAD_TYPE: &'static str = "edge_traction";

    /// Returns the loaded element ID.
    #[must_use]
    pub fn element_id(&self) -> usize {
        self.element_id
    }

    /// Returns the two node IDs identifying the loaded edge.
    #[must_use]
    pub fn edge_node_ids(&self) -> [usize; 2] {
        self.edge_node_ids
    }

    /// Returns the coordinate system used by the traction components.
    #[must_use]
    pub fn coordinate_system(&self) -> LoadCoordinateSystem2D {
        self.coordinate_system
    }

    /// Returns the x component of the traction.
    #[must_use]
    pub fn x_component(&self) -> f64 {
        self.x_component
    }

    /// Returns the y component of the traction.
    #[must_use]
    pub fn y_component(&self) -> f64 {
        self.y_component
    }

    pub(crate) fn global_components(&self, cosine: f64, sine: f64) -> (f64, f64) {
        match self.coordinate_system {
            LoadCoordinateSystem2D::Global => (self.x_component, self.y_component),
            LoadCoordinateSystem2D::Local => (
                cosine * self.x_component - sine * self.y_component,
                sine * self.x_component + cosine * self.y_component,
            ),
        }
    }
}

/// Uniform body force applied over a plane-stress T3 element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BodyForce2D {
    element_id: usize,
    x_component: f64,
    y_component: f64,
}

impl BodyForce2D {
    /// Creates a uniform body force.
    ///
    /// The components are force per unit volume in the model's global x/y
    /// coordinate system.
    pub fn new(element_id: usize, x_component: f64, y_component: f64) -> Result<Self, FemError> {
        if !x_component.is_finite() {
            return Err(FemError::InvalidElementLoadValue {
                element_id,
                load_type: Self::LOAD_TYPE,
                component: "x component",
                value: x_component,
            });
        }

        if !y_component.is_finite() {
            return Err(FemError::InvalidElementLoadValue {
                element_id,
                load_type: Self::LOAD_TYPE,
                component: "y component",
                value: y_component,
            });
        }

        Ok(Self { element_id, x_component, y_component })
    }

    pub(crate) const LOAD_TYPE: &'static str = "body_force";

    /// Returns the loaded element ID.
    #[must_use]
    pub fn element_id(&self) -> usize {
        self.element_id
    }

    /// Returns the global x component of the body force.
    #[must_use]
    pub fn x_component(&self) -> f64 {
        self.x_component
    }

    /// Returns the global y component of the body force.
    #[must_use]
    pub fn y_component(&self) -> f64 {
        self.y_component
    }
}

/// Self-weight load applied over a plane-stress T3 element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SelfWeight2D {
    element_id: usize,
    x_acceleration: f64,
    y_acceleration: f64,
}

impl SelfWeight2D {
    /// Creates a self-weight load.
    ///
    /// The acceleration components are defined in the model's global x/y
    /// coordinate system. During assembly they are multiplied by the loaded
    /// element material density.
    pub fn new(element_id: usize, x_acceleration: f64, y_acceleration: f64) -> Result<Self, FemError> {
        if !x_acceleration.is_finite() {
            return Err(FemError::InvalidElementLoadValue {
                element_id,
                load_type: Self::LOAD_TYPE,
                component: "x acceleration",
                value: x_acceleration,
            });
        }

        if !y_acceleration.is_finite() {
            return Err(FemError::InvalidElementLoadValue {
                element_id,
                load_type: Self::LOAD_TYPE,
                component: "y acceleration",
                value: y_acceleration,
            });
        }

        Ok(Self { element_id, x_acceleration, y_acceleration })
    }

    pub(crate) const LOAD_TYPE: &'static str = "self_weight";

    /// Returns the loaded element ID.
    #[must_use]
    pub fn element_id(&self) -> usize {
        self.element_id
    }

    /// Returns the global x acceleration component.
    #[must_use]
    pub fn x_acceleration(&self) -> f64 {
        self.x_acceleration
    }

    /// Returns the global y acceleration component.
    #[must_use]
    pub fn y_acceleration(&self) -> f64 {
        self.y_acceleration
    }
}

/// Loads that act over an element rather than directly at a node.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ElementLoad2D {
    /// Uniform line load along a beam element.
    BeamUniformLine(BeamUniformLineLoad2D),

    /// Uniform traction along one T3 element edge.
    EdgeTraction(EdgeTraction2D),

    /// Uniform body force over one T3 element.
    BodyForce(BodyForce2D),

    /// Self-weight load over one T3 element.
    SelfWeight(SelfWeight2D),
}

impl ElementLoad2D {
    /// Returns the loaded element ID.
    #[must_use]
    pub fn element_id(&self) -> usize {
        match self {
            Self::BeamUniformLine(load) => load.element_id(),
            Self::EdgeTraction(load) => load.element_id(),
            Self::BodyForce(load) => load.element_id(),
            Self::SelfWeight(load) => load.element_id(),
        }
    }

    /// Returns a stable load type name for diagnostics.
    #[must_use]
    pub fn load_type(&self) -> &'static str {
        match self {
            Self::BeamUniformLine(_) => BeamUniformLineLoad2D::LOAD_TYPE,
            Self::EdgeTraction(_) => EdgeTraction2D::LOAD_TYPE,
            Self::BodyForce(_) => BodyForce2D::LOAD_TYPE,
            Self::SelfWeight(_) => SelfWeight2D::LOAD_TYPE,
        }
    }

    /// Returns the element type expected by this load.
    #[must_use]
    pub fn expected_element_type(&self) -> &'static str {
        match self {
            Self::BeamUniformLine(_) => "beam",
            Self::EdgeTraction(_) => "triangle_t3",
            Self::BodyForce(_) => "triangle_t3",
            Self::SelfWeight(_) => "triangle_t3",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{
        BeamUniformLineLoad2D, BodyForce2D, EdgeTraction2D, ElementLoad2D, LoadCoordinateSystem2D, SelfWeight2D,
    };
    use crate::FemError;

    #[test]
    fn creates_uniform_beam_line_load_with_valid_data() {
        let load = BeamUniformLineLoad2D::new(10, LoadCoordinateSystem2D::Local, 2.0, -5.0)
            .expect("valid load should be created");

        assert_eq!(load.element_id(), 10);
        assert_eq!(load.coordinate_system(), LoadCoordinateSystem2D::Local);
        assert_eq!(load.x_component(), 2.0);
        assert_eq!(load.y_component(), -5.0);
    }

    #[test]
    fn rejects_non_finite_uniform_beam_line_load_values() {
        let result = BeamUniformLineLoad2D::new(10, LoadCoordinateSystem2D::Local, f64::NAN, 0.0);

        assert!(matches!(
            result,
            Err(FemError::InvalidElementLoadValue {
                element_id: 10,
                load_type: "beam_uniform_line",
                component: "x component",
                value
            }) if value.is_nan()
        ));

        let result = BeamUniformLineLoad2D::new(20, LoadCoordinateSystem2D::Global, 0.0, f64::INFINITY);

        assert!(matches!(
            result,
            Err(FemError::InvalidElementLoadValue {
                element_id: 20,
                load_type: "beam_uniform_line",
                component: "y component",
                value
            }) if value.is_infinite()
        ));
    }

    #[test]
    fn converts_global_components_to_local_components() {
        let load = BeamUniformLineLoad2D::new(10, LoadCoordinateSystem2D::Global, 0.0, -10.0)
            .expect("valid load should be created");
        let cosine = 0.6;
        let sine = 0.8;

        let (x_component, y_component) = load.local_components(cosine, sine);

        assert!((x_component + 8.0).abs() < 1e-12);
        assert!((y_component + 6.0).abs() < 1e-12);
    }

    #[test]
    fn calculates_local_equivalent_nodal_load() {
        let load = BeamUniformLineLoad2D::new(10, LoadCoordinateSystem2D::Local, 2.0, -3.0)
            .expect("valid load should be created");

        let actual = load.local_equivalent_nodal_load(4.0, 1.0, 0.0);

        assert_eq!(actual, [4.0, -6.0, -4.0, 4.0, -6.0, 4.0]);
    }

    #[test]
    fn element_load_reports_diagnostics_data() {
        let load = BeamUniformLineLoad2D::new(10, LoadCoordinateSystem2D::Local, 0.0, -1.0)
            .expect("valid load should be created");
        let load = ElementLoad2D::BeamUniformLine(load);

        assert_eq!(load.element_id(), 10);
        assert_eq!(load.load_type(), "beam_uniform_line");
        assert_eq!(load.expected_element_type(), "beam");
    }

    #[test]
    fn creates_edge_traction_with_valid_data() {
        let load = EdgeTraction2D::new(10, [2, 3], LoadCoordinateSystem2D::Global, 5.0, -2.0)
            .expect("valid load should be created");

        assert_eq!(load.element_id(), 10);
        assert_eq!(load.edge_node_ids(), [2, 3]);
        assert_eq!(load.coordinate_system(), LoadCoordinateSystem2D::Global);
        assert_eq!(load.x_component(), 5.0);
        assert_eq!(load.y_component(), -2.0);
    }

    #[test]
    fn rejects_edge_traction_with_repeated_edge_node_ids() {
        let result = EdgeTraction2D::new(10, [2, 2], LoadCoordinateSystem2D::Global, 5.0, -2.0);

        assert!(matches!(
            result,
            Err(FemError::InvalidElementLoadEdge {
                element_id: 10,
                load_type: "edge_traction",
                node_ids,
                expected: "two distinct edge node IDs",
            }) if node_ids == vec![2, 2]
        ));
    }

    #[test]
    fn rejects_non_finite_edge_traction_values() {
        let result = EdgeTraction2D::new(10, [2, 3], LoadCoordinateSystem2D::Global, f64::NAN, 0.0);

        assert!(matches!(
            result,
            Err(FemError::InvalidElementLoadValue {
                element_id: 10,
                load_type: "edge_traction",
                component: "x component",
                value
            }) if value.is_nan()
        ));

        let result = EdgeTraction2D::new(20, [2, 3], LoadCoordinateSystem2D::Local, 0.0, f64::INFINITY);

        assert!(matches!(
            result,
            Err(FemError::InvalidElementLoadValue {
                element_id: 20,
                load_type: "edge_traction",
                component: "y component",
                value
            }) if value.is_infinite()
        ));
    }

    #[test]
    fn converts_local_edge_traction_components_to_global_components() {
        let load = EdgeTraction2D::new(10, [2, 3], LoadCoordinateSystem2D::Local, 10.0, -5.0)
            .expect("valid load should be created");
        let cosine = 0.6;
        let sine = 0.8;

        let (x_component, y_component) = load.global_components(cosine, sine);

        assert!((x_component - 10.0).abs() < 1e-12);
        assert!((y_component - 5.0).abs() < 1e-12);
    }

    #[test]
    fn edge_traction_reports_diagnostics_data() {
        let load = EdgeTraction2D::new(10, [2, 3], LoadCoordinateSystem2D::Global, 0.0, -1.0)
            .expect("valid load should be created");
        let load = ElementLoad2D::EdgeTraction(load);

        assert_eq!(load.element_id(), 10);
        assert_eq!(load.load_type(), "edge_traction");
        assert_eq!(load.expected_element_type(), "triangle_t3");
    }

    #[test]
    fn creates_body_force_with_valid_data() {
        let load = BodyForce2D::new(10, 5.0, -2.0).expect("valid load should be created");

        assert_eq!(load.element_id(), 10);
        assert_eq!(load.x_component(), 5.0);
        assert_eq!(load.y_component(), -2.0);
    }

    #[test]
    fn rejects_non_finite_body_force_values() {
        let result = BodyForce2D::new(10, f64::NAN, 0.0);

        assert!(matches!(
            result,
            Err(FemError::InvalidElementLoadValue {
                element_id: 10,
                load_type: "body_force",
                component: "x component",
                value
            }) if value.is_nan()
        ));

        let result = BodyForce2D::new(20, 0.0, f64::INFINITY);

        assert!(matches!(
            result,
            Err(FemError::InvalidElementLoadValue {
                element_id: 20,
                load_type: "body_force",
                component: "y component",
                value
            }) if value.is_infinite()
        ));
    }

    #[test]
    fn body_force_reports_diagnostics_data() {
        let load = BodyForce2D::new(10, 0.0, -1.0).expect("valid load should be created");
        let load = ElementLoad2D::BodyForce(load);

        assert_eq!(load.element_id(), 10);
        assert_eq!(load.load_type(), "body_force");
        assert_eq!(load.expected_element_type(), "triangle_t3");
    }

    #[test]
    fn creates_self_weight_with_valid_data() {
        let load = SelfWeight2D::new(10, 0.0, -9.81).expect("valid load should be created");

        assert_eq!(load.element_id(), 10);
        assert_eq!(load.x_acceleration(), 0.0);
        assert_eq!(load.y_acceleration(), -9.81);
    }

    #[test]
    fn rejects_non_finite_self_weight_accelerations() {
        let result = SelfWeight2D::new(10, f64::NAN, 0.0);

        assert!(matches!(
            result,
            Err(FemError::InvalidElementLoadValue {
                element_id: 10,
                load_type: "self_weight",
                component: "x acceleration",
                value
            }) if value.is_nan()
        ));

        let result = SelfWeight2D::new(20, 0.0, f64::INFINITY);

        assert!(matches!(
            result,
            Err(FemError::InvalidElementLoadValue {
                element_id: 20,
                load_type: "self_weight",
                component: "y acceleration",
                value
            }) if value.is_infinite()
        ));
    }

    #[test]
    fn self_weight_reports_diagnostics_data() {
        let load = SelfWeight2D::new(10, 0.0, -9.81).expect("valid load should be created");
        let load = ElementLoad2D::SelfWeight(load);

        assert_eq!(load.element_id(), 10);
        assert_eq!(load.load_type(), "self_weight");
        assert_eq!(load.expected_element_type(), "triangle_t3");
    }
}
