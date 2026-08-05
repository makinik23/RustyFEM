//! Material data for a two-dimensional finite element model.

use super::material_model::Material2D;
use crate::FemError;

/// Material ID used by the compatibility single-material API.
pub const DEFAULT_MATERIAL_ID: usize = 0;

/// Stores materials used by a 2D model.
#[derive(Default)]
pub struct Materials2D {
    materials: Vec<(usize, Material2D)>,
}

impl Materials2D {
    /// Creates an empty material collection.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a material with a model-unique ID.
    pub fn add_material(&mut self, material_id: usize, material: Material2D) -> Result<(), FemError> {
        if self.contains_material_id(material_id) {
            return Err(FemError::DuplicateId { entity: "material", id: material_id });
        }

        self.materials.push((material_id, material));

        Ok(())
    }

    /// Sets the default material used by the current single-material model API.
    pub(crate) fn set_default_material(&mut self, material: Material2D) {
        match self.materials.iter_mut().find(|(id, _)| *id == DEFAULT_MATERIAL_ID) {
            Some((_, existing)) => *existing = material,
            None => self.materials.push((DEFAULT_MATERIAL_ID, material)),
        }
    }

    /// Returns the default material, if one was set.
    #[must_use]
    pub fn default_material(&self) -> Option<&Material2D> {
        self.materials.iter().find(|(id, _)| *id == DEFAULT_MATERIAL_ID).map(|(_, material)| material)
    }

    /// Returns a material by ID.
    pub fn material(&self, material_id: usize) -> Result<&Material2D, FemError> {
        self.materials
            .iter()
            .find(|(id, _)| *id == material_id)
            .map(|(_, material)| material)
            .ok_or(FemError::UnknownId { entity: "material", id: material_id })
    }

    /// Returns all materials in insertion order.
    #[must_use]
    pub fn materials(&self) -> &[(usize, Material2D)] {
        &self.materials
    }

    /// Checks whether a material ID already exists.
    #[must_use]
    pub(crate) fn contains_material_id(&self, material_id: usize) -> bool {
        self.materials.iter().any(|(id, _)| *id == material_id)
    }
}

#[cfg(test)]
mod tests {
    use super::{DEFAULT_MATERIAL_ID, Materials2D};
    use crate::FemError;
    use crate::model::Material2D;

    #[test]
    fn creates_empty_materials() {
        let materials = Materials2D::new();

        assert!(materials.default_material().is_none());
        assert!(materials.materials().is_empty());
    }

    #[test]
    fn stores_default_material() {
        let mut materials = Materials2D::new();
        let material = Material2D::new(210e9, 0.3, 7800.0).expect("valid material");

        materials.set_default_material(material);

        assert_eq!(materials.default_material(), Some(&material));
        assert_eq!(materials.material(DEFAULT_MATERIAL_ID).expect("default material should exist"), &material);
    }

    #[test]
    fn replaces_default_material() {
        let mut materials = Materials2D::new();
        let first = Material2D::new(210e9, 0.3, 7800.0).expect("valid material");
        let second = Material2D::new(70e9, 0.33, 2700.0).expect("valid material");

        materials.set_default_material(first);
        materials.set_default_material(second);

        assert_eq!(materials.default_material(), Some(&second));
    }

    #[test]
    fn stores_materials_in_insertion_order() {
        let mut materials = Materials2D::new();
        let steel = Material2D::new(210e9, 0.3, 7800.0).expect("valid material");
        let aluminium = Material2D::new(70e9, 0.33, 2700.0).expect("valid material");

        materials.add_material(10, steel).expect("steel should be added");
        materials.add_material(20, aluminium).expect("aluminium should be added");

        assert_eq!(materials.materials(), &[(10, steel), (20, aluminium)]);
        assert_eq!(materials.material(10).expect("steel should exist"), &steel);
        assert_eq!(materials.material(20).expect("aluminium should exist"), &aluminium);
    }

    #[test]
    fn rejects_duplicate_material_ids() {
        let mut materials = Materials2D::new();
        let material = Material2D::new(210e9, 0.3, 7800.0).expect("valid material");

        materials.add_material(10, material).expect("first material should be added");
        let result = materials.add_material(10, material);

        assert!(matches!(result, Err(FemError::DuplicateId { entity: "material", id: 10 })));
    }

    #[test]
    fn rejects_unknown_material_ids() {
        let materials = Materials2D::new();

        assert!(matches!(materials.material(99), Err(FemError::UnknownId { entity: "material", id: 99 })));
    }
}
