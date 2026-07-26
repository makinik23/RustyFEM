//! Defines a single node by it's spacial coordinates and ID.

use crate::FemError;

/// 2D node, defined by its ID and coordinates (x, y).
pub struct Node2D {
    pub id: usize,
    pub x: f64,
    pub y: f64,
}

impl Node2D {
    /// Creates a new 2D node with the given ID and coordinates.
    pub fn new(id: usize, x: f64, y: f64) -> Result<Self, FemError> {
        if !x.is_finite() {
            return Err(FemError::InvalidNodeCoordinate { node_id: id, axis: "x", value: x });
        }

        if !y.is_finite() {
            return Err(FemError::InvalidNodeCoordinate { node_id: id, axis: "y", value: y });
        }

        Ok(Self { id, x, y })
    }

    /// Returns the ID of the node.
    #[must_use]
    pub fn id(&self) -> usize {
        self.id
    }

    /// Returns the x-coordinate of the node.
    #[must_use]
    pub fn x(&self) -> f64 {
        self.x
    }

    /// Returns the y-coordinate of the node.
    #[must_use]
    pub fn y(&self) -> f64 {
        self.y
    }

    /// Returns the coordinates of the node as a tuple (x, y).
    #[must_use]
    pub fn coordinates(&self) -> (f64, f64) {
        (self.x, self.y)
    }
}

#[cfg(test)]
mod tests {
    use super::Node2D;
    use crate::FemError;

    #[test]
    fn creates_node_with_valid_data() {
        let node = Node2D::new(7, 1.5, -2.0).expect("valid node should be created");

        assert_eq!(node.id(), 7);
        assert_eq!(node.x(), 1.5);
        assert_eq!(node.y(), -2.0);
        assert_eq!(node.coordinates(), (1.5, -2.0));
    }

    #[test]
    fn rejects_nan_x_coordinate() {
        let result = Node2D::new(1, f64::NAN, 0.0);

        assert!(matches!(
            result,
            Err(FemError::InvalidNodeCoordinate
            {
                node_id: 1,
                axis: "x",
                value,
            }) if value.is_nan()
        ));
    }

    #[test]
    fn rejects_nan_y_coordinate() {
        let result = Node2D::new(1, 0.0, f64::NAN);

        assert!(matches!(
            result,
            Err(FemError::InvalidNodeCoordinate
            {
                node_id: 1,
                axis: "y",
                value,
            }) if value.is_nan()
        ));
    }
}
