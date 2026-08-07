use std::io::{self, Write};

use clap::{Parser, ValueEnum};
use rusty_fem::FemError;
use rusty_fem::analysis::iterative_solver::CgTerminationReason;
use rusty_fem::analysis::solver::{AnalysisResult2D, solve_with_settings};
use rusty_fem::analysis::{ElementResponse2D, recover_beam_section_response, recover_model_responses};
use rusty_fem::elements::{Beam2D, Element2D, TriangleT3, Truss2D};
use rusty_fem::model::{
    AnalysisSpace, BeamSection2D, BeamUniformLineLoad2D, BodyForce2D, DisplacementConstraint2D, Dof2D, DofNumbering2D,
    EdgeTraction2D, ElementLoad2D, LoadCoordinateSystem2D, Material2D, Model2D, NodalLoad2D, Node2D,
    PlaneStressSection2D, Section2D, SelfWeight2D, SolverKind2D, TrussSection2D,
};

#[derive(Debug, Parser)]
#[command(name = "rusty-fem", about = "Educational finite element method solver")]
struct Cli {
    #[arg(long, value_name = "2D|3D")]
    space: Option<AnalysisSpace>,

    /// Selects the linear solver used for the analysis.
    #[arg(long, value_enum, default_value_t = SolverKind::Dense)]
    solver: SolverKind,

    /// Relative residual tolerance used by the sparse iterative solver.
    #[arg(long, default_value_t = 1e-10, value_name = "TOLERANCE")]
    cg_tolerance: f64,

    /// Maximum number of CG iterations used by the sparse solver.
    #[arg(long, default_value_t = 1_000, value_name = "ITERATIONS")]
    cg_max_iterations: usize,
}

/// Selects between the reference dense solver and the sparse iterative solver.
#[derive(Debug, Clone, Copy, ValueEnum)]
enum SolverKind {
    /// Solves the system with a dense LU decomposition.
    Dense,

    /// Solves the system with CSR storage and Conjugate Gradient.
    Sparse,
}

impl From<SolverKind> for SolverKind2D {
    fn from(value: SolverKind) -> Self {
        match value {
            SolverKind::Dense => Self::Dense,
            SolverKind::Sparse => Self::Sparse,
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    println!("RustyFEM interactive session");

    let space = match cli.space {
        Some(space) => space,
        None => read_analysis_space()?,
    };

    println!("Selected analysis space: {:?}", space);
    println!("Spatial dimension: {}", space.spatial_dimension());

    if space == AnalysisSpace::ThreeDimensional {
        println!("3D logic is not implemented yet.");

        return Ok(());
    }

    let mut model = Model2D::new();
    model.analysis_settings_mut().set_solver(cli.solver.into());
    model.analysis_settings_mut().set_cg_tolerance(cli.cg_tolerance)?;
    model.analysis_settings_mut().set_cg_max_iterations(cli.cg_max_iterations)?;

    read_materials(&mut model)?;
    read_sections(&mut model)?;
    read_nodes(&mut model)?;
    read_constraints(&mut model)?;
    read_elements(&mut model)?;
    read_loads(&mut model)?;

    println!();
    println!("Model summary:");
    println!("  materials: {}", model.materials().materials().len());
    println!("  sections: {}", model.sections().sections().len());
    println!("  nodes: {}", model.nodes().len());
    println!("  constraints: {}", model.constraints().len());
    println!("  elements: {}", model.elements().len());
    println!(
        "  loads: {} ({} nodal, {} element)",
        model.loads().len() + model.element_loads().len(),
        model.loads().len(),
        model.element_loads().len()
    );

    println!();
    let result = solve_with_settings(&model);

    match result {
        Ok(result) => print_analysis_results(&model, &result)?,
        Err(error) => println!("Could not solve model: {error}"),
    }

    Ok(())
}

fn prompt_line(prompt: &str) -> io::Result<Option<String>> {
    print!("{prompt}");
    io::stdout().flush()?;

    let mut line = String::new();
    let bytes_read = io::stdin().read_line(&mut line)?;

    if bytes_read == 0 {
        return Ok(None);
    }

    Ok(Some(line.trim().to_owned()))
}

fn read_analysis_space() -> io::Result<AnalysisSpace> {
    loop {
        let Some(line) = prompt_line("Analysis space [2d/3d]: ")? else {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "input ended"));
        };

        match line.parse::<AnalysisSpace>() {
            Ok(space) => return Ok(space),
            Err(error) => println!("Invalid analysis space: {error}"),
        }
    }
}

fn read_materials(model: &mut Model2D) -> io::Result<()> {
    println!();
    println!("Enter materials as: ID YOUNG_MODULUS POISSON_RATIO DENSITY");
    println!("Type 'done' when finished.");

    loop {
        let Some(line) = prompt_line("material> ")? else {
            return Ok(());
        };

        if line.eq_ignore_ascii_case("done") {
            return Ok(());
        }

        match parse_material_line(&line) {
            Ok((material_id, material)) => match model.add_material(material_id, material) {
                Ok(()) => println!("Material {material_id} added."),
                Err(error) => println!("Could not add material: {error}"),
            },
            Err(error) => println!("Invalid material: {error}"),
        }
    }
}

fn read_sections(model: &mut Model2D) -> io::Result<()> {
    println!();
    println!("Enter sections as:");
    println!("  truss ID AREA");
    println!("  beam ID AREA I HEIGHT");
    println!("  plane_stress ID THICKNESS");
    println!("Type 'done' when finished.");

    loop {
        let Some(line) = prompt_line("section> ")? else {
            return Ok(());
        };

        if line.eq_ignore_ascii_case("done") {
            return Ok(());
        }

        match parse_section_line(&line) {
            Ok((section_id, section)) => match model.add_section(section_id, section) {
                Ok(()) => println!("Section {section_id} added."),
                Err(error) => println!("Could not add section: {error}"),
            },
            Err(error) => println!("Invalid section: {error}"),
        }
    }
}

fn read_nodes(model: &mut Model2D) -> io::Result<()> {
    println!();
    println!("Enter nodes as: ID X Y");
    println!("Type 'done' when finished.");

    loop {
        let Some(line) = prompt_line("node> ")? else {
            return Ok(());
        };

        if line.eq_ignore_ascii_case("done") {
            return Ok(());
        }

        match parse_node_line(&line) {
            Ok(node) => match model.add_node(node) {
                Ok(()) => println!("Node added."),
                Err(error) => println!("Could not add node: {error}"),
            },
            Err(error) => println!("Invalid node: {error}"),
        }
    }
}

fn read_constraints(model: &mut Model2D) -> io::Result<()> {
    println!();
    println!("Enter constraints as: NODE_ID DOF VALUE");
    println!("DOF must be Ux, Uy, or Rz.");
    println!("Type 'done' when finished.");

    loop {
        let Some(line) = prompt_line("constraint> ")? else {
            return Ok(());
        };

        if line.eq_ignore_ascii_case("done") {
            return Ok(());
        }

        match parse_constraint_line(&line) {
            Ok(constraint) => match model.add_constraint(constraint) {
                Ok(()) => println!("Constraint added."),
                Err(error) => println!("Could not add constraint: {error}"),
            },
            Err(error) => println!("Invalid constraint: {error}"),
        }
    }
}

fn read_elements(model: &mut Model2D) -> io::Result<()> {
    println!();
    println!("Enter elements as:");
    println!("  truss ID NODE_1 NODE_2 MATERIAL_ID SECTION_ID");
    println!("  beam ID NODE_1 NODE_2 MATERIAL_ID SECTION_ID");
    println!("  triangle ID NODE_1 NODE_2 NODE_3 MATERIAL_ID SECTION_ID");
    println!("Type 'done' when finished.");

    loop {
        let Some(line) = prompt_line("element> ")? else {
            return Ok(());
        };

        if line.eq_ignore_ascii_case("done") {
            return Ok(());
        }

        match parse_element_line(&line) {
            Ok(element) => {
                let element_id = element.id();
                let interpolation = element.interpolation();
                let dof_count = element.dof_count();

                match model.add_element(element) {
                    Ok(()) => println!("Element {element_id} added ({interpolation:?}, {dof_count} DOFs)."),
                    Err(error) => println!("Could not add element: {error}"),
                }
            }
            Err(error) => println!("Invalid element: {error}"),
        }
    }
}

fn read_loads(model: &mut Model2D) -> io::Result<()> {
    println!();
    println!("Enter nodal loads as: NODE_ID DOF VALUE");
    println!("Enter uniform beam line loads as: beam_uniform ELEMENT_ID local|global QX QY");
    println!("Enter T3 edge tractions as: edge_traction ELEMENT_ID NODE_A NODE_B local|global TX TY");
    println!("Enter T3 body forces as: body_force ELEMENT_ID global BX BY");
    println!("Enter T3 self-weight loads as: self_weight ELEMENT_ID global AX AY");
    println!(
        "DOF must be Ux, Uy, or Rz. Beam loads are force per unit length; T3 loads use force per unit area or volume."
    );
    println!("Type 'done' when finished.");

    loop {
        let Some(line) = prompt_line("load> ")? else {
            return Ok(());
        };

        if line.eq_ignore_ascii_case("done") {
            return Ok(());
        }

        match parse_load_line(&line) {
            Ok(LoadInput2D::Nodal(load)) => match model.add_load(load) {
                Ok(()) => println!("Nodal load added."),
                Err(error) => println!("Could not add load: {error}"),
            },
            Ok(LoadInput2D::Element(load)) => match model.add_element_load(load) {
                Ok(()) => println!("Element load added."),
                Err(error) => println!("Could not add load: {error}"),
            },
            Err(error) => println!("Invalid load: {error}"),
        }
    }
}

fn print_analysis_results(model: &Model2D, result: &AnalysisResult2D) -> Result<(), Box<dyn std::error::Error>> {
    let numbering = DofNumbering2D::from_model(model)?;
    let ordered_dofs = [Dof2D::Ux, Dof2D::Uy, Dof2D::Rz];

    if let Some(report) = result.solver_report() {
        let status = match report.termination_reason {
            CgTerminationReason::Converged => "converged",
            CgTerminationReason::MaxIterations => "maximum iterations reached",
            CgTerminationReason::Stagnated => "stagnated",
        };

        println!("Solver report:");
        println!("  iterations: {}", report.iterations);
        println!("  residual norm: {:.6e}", report.residual_norm);
        println!("  relative residual norm: {:.6e}", report.relative_residual_norm);
        println!("  status: {status}");
    }

    println!("Displacements:");

    for node in model.nodes() {
        for dof in ordered_dofs {
            if let Ok(index) = numbering.index(node.id(), dof) {
                println!("  node {} {} = {:.12}", node.id(), dof.name(), result.displacements()[index]);
            }
        }
    }

    println!("Reactions:");

    for constraint in model.constraints() {
        let index = numbering.index(constraint.node_id(), constraint.dof())?;

        println!("  node {} {} = {:.12}", constraint.node_id(), constraint.dof().name(), result.reactions()[index]);
    }

    println!("Element responses:");

    for (element_id, response) in recover_model_responses(model, result.displacements())? {
        print_element_response(element_id, response);
    }

    print_beam_section_results(model, result)?;

    Ok(())
}

fn print_element_response(element_id: usize, response: ElementResponse2D) {
    match response {
        ElementResponse2D::Truss(response) => {
            println!(
                "  element {element_id} truss: strain = {:.12}, stress = {:.12}, axial_force = {:.12}",
                response.strain(),
                response.stress(),
                response.axial_force(),
            );
        }
        ElementResponse2D::Beam(response) => {
            println!(
                "  element {element_id} beam: [N1, V1, M1, N2, V2, M2] = [{:.12}, {:.12}, {:.12}, {:.12}, {:.12}, {:.12}]",
                response.first_axial_force(),
                response.first_shear_force(),
                response.first_bending_moment(),
                response.second_axial_force(),
                response.second_shear_force(),
                response.second_bending_moment(),
            );
        }
        ElementResponse2D::Triangle(response) => {
            let strain = response.strain();
            let stress = response.stress();

            println!(
                "  element {element_id} triangle: strain = [{:.12}, {:.12}, {:.12}], stress = [{:.12}, {:.12}, {:.12}], von_mises = {:.12}",
                strain[0],
                strain[1],
                strain[2],
                stress[0],
                stress[1],
                stress[2],
                response.von_mises_stress(),
            );
        }
    }
}

fn print_beam_section_results(model: &Model2D, result: &AnalysisResult2D) -> Result<(), Box<dyn std::error::Error>> {
    for element in model.elements() {
        let Element2D::Beam(beam) = element else {
            continue;
        };

        let node_ids = element.node_ids();
        let first_node = model
            .nodes()
            .iter()
            .find(|node| node.id() == node_ids[0])
            .ok_or(FemError::UnknownId { entity: "node", id: node_ids[0] })?;
        let second_node = model
            .nodes()
            .iter()
            .find(|node| node.id() == node_ids[1])
            .ok_or(FemError::UnknownId { entity: "node", id: node_ids[1] })?;
        let length = beam.length(first_node, second_node)?;

        println!("Beam {} section results:", element.id());

        for position in [0.0, length / 2.0, length] {
            let response = recover_beam_section_response(model, beam, result.displacements(), position)?;

            let section = model.beam_section(beam.section_id())?;

            match section.section_height() {
                Some(height) => println!(
                    "  x = {:.12}: curvature = {:.12}, bending_moment = {:.12}, sigma(y=+h/2) = {:.12}, sigma(y=-h/2) = {:.12}",
                    position,
                    response.curvature(),
                    response.bending_moment(),
                    response.normal_stress(height / 2.0),
                    response.normal_stress(-height / 2.0),
                ),
                None => println!(
                    "  x = {:.12}: curvature = {:.12}, bending_moment = {:.12}, fiber_stress = unavailable (section height not provided)",
                    position,
                    response.curvature(),
                    response.bending_moment(),
                ),
            }
        }
    }

    Ok(())
}

fn parse_node_line(line: &str) -> Result<Node2D, String> {
    let parts: Vec<&str> = line.split_whitespace().collect();

    if parts.len() != 3 {
        return Err("expected: ID X Y".to_owned());
    }

    let id = parse_usize(parts[0], "node ID")?;
    let x = parse_f64(parts[1], "x-coordinate")?;
    let y = parse_f64(parts[2], "y-coordinate")?;

    Node2D::new(id, x, y).map_err(|error| error.to_string())
}

fn parse_constraint_line(line: &str) -> Result<DisplacementConstraint2D, String> {
    let parts: Vec<&str> = line.split_whitespace().collect();

    if parts.len() != 3 {
        return Err("expected: NODE_ID DOF VALUE".to_owned());
    }

    let node_id = parse_usize(parts[0], "node ID")?;
    let dof = parse_dof(parts[1])?;
    let displacement = parse_f64(parts[2], "displacement")?;

    DisplacementConstraint2D::new(node_id, dof, displacement).map_err(|error| error.to_string())
}

fn parse_material_line(line: &str) -> Result<(usize, Material2D), String> {
    let parts: Vec<&str> = line.split_whitespace().collect();

    if parts.len() != 4 {
        return Err("expected: ID YOUNG_MODULUS POISSON_RATIO DENSITY".to_owned());
    }

    let id = parse_usize(parts[0], "material ID")?;
    let young_modulus = parse_f64(parts[1], "Young's modulus")?;
    let poisson_ratio = parse_f64(parts[2], "Poisson's ratio")?;
    let density = parse_f64(parts[3], "density")?;

    let material = Material2D::new(young_modulus, poisson_ratio, density).map_err(|error| error.to_string())?;

    Ok((id, material))
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum LoadInput2D {
    Nodal(NodalLoad2D),
    Element(ElementLoad2D),
}

fn parse_load_line(line: &str) -> Result<LoadInput2D, String> {
    let parts: Vec<&str> = line.split_whitespace().collect();

    if parts.is_empty() {
        return Err("expected a load".to_owned());
    }

    match parts[0].to_ascii_lowercase().as_str() {
        "beam_uniform" | "beam_uniform_line" => parse_beam_uniform_line_load(&parts).map(LoadInput2D::Element),
        "edge_traction" => parse_edge_traction_load(&parts).map(LoadInput2D::Element),
        "body_force" => parse_body_force_load(&parts).map(LoadInput2D::Element),
        "self_weight" => parse_self_weight_load(&parts).map(LoadInput2D::Element),
        "nodal" => parse_nodal_load_parts(&parts[1..]).map(LoadInput2D::Nodal),
        _ => parse_nodal_load_parts(&parts).map(LoadInput2D::Nodal),
    }
}

fn parse_nodal_load_parts(parts: &[&str]) -> Result<NodalLoad2D, String> {
    if parts.len() != 3 {
        return Err("expected: NODE_ID DOF VALUE".to_owned());
    }

    let node_id = parse_usize(parts[0], "node ID")?;
    let dof = parse_dof(parts[1])?;
    let value = parse_f64(parts[2], "load value")?;

    NodalLoad2D::new(node_id, dof, value).map_err(|error| error.to_string())
}

fn parse_beam_uniform_line_load(parts: &[&str]) -> Result<ElementLoad2D, String> {
    if parts.len() != 5 {
        return Err("expected: beam_uniform ELEMENT_ID local|global QX QY".to_owned());
    }

    let element_id = parse_usize(parts[1], "element ID")?;
    let coordinate_system = parse_load_coordinate_system(parts[2])?;
    let x_component = parse_f64(parts[3], "line load x component")?;
    let y_component = parse_f64(parts[4], "line load y component")?;
    let load = BeamUniformLineLoad2D::new(element_id, coordinate_system, x_component, y_component)
        .map_err(|error| error.to_string())?;

    Ok(ElementLoad2D::BeamUniformLine(load))
}

fn parse_edge_traction_load(parts: &[&str]) -> Result<ElementLoad2D, String> {
    if parts.len() != 7 {
        return Err("expected: edge_traction ELEMENT_ID NODE_A NODE_B local|global TX TY".to_owned());
    }

    let element_id = parse_usize(parts[1], "element ID")?;
    let first_node_id = parse_usize(parts[2], "first edge node ID")?;
    let second_node_id = parse_usize(parts[3], "second edge node ID")?;
    let coordinate_system = parse_load_coordinate_system(parts[4])?;
    let x_component = parse_f64(parts[5], "edge traction x component")?;
    let y_component = parse_f64(parts[6], "edge traction y component")?;
    let load =
        EdgeTraction2D::new(element_id, [first_node_id, second_node_id], coordinate_system, x_component, y_component)
            .map_err(|error| error.to_string())?;

    Ok(ElementLoad2D::EdgeTraction(load))
}

fn parse_body_force_load(parts: &[&str]) -> Result<ElementLoad2D, String> {
    let (element_id, x_component, y_component) = match parts {
        [_, element_id, coordinate_system, x_component, y_component] => {
            if !coordinate_system.eq_ignore_ascii_case("global") {
                return Err(format!("unknown body force coordinate system: {coordinate_system}"));
            }

            (
                parse_usize(element_id, "element ID")?,
                parse_f64(x_component, "body force x component")?,
                parse_f64(y_component, "body force y component")?,
            )
        }
        [_, element_id, x_component, y_component] => (
            parse_usize(element_id, "element ID")?,
            parse_f64(x_component, "body force x component")?,
            parse_f64(y_component, "body force y component")?,
        ),
        _ => return Err("expected: body_force ELEMENT_ID global BX BY".to_owned()),
    };
    let load = BodyForce2D::new(element_id, x_component, y_component).map_err(|error| error.to_string())?;

    Ok(ElementLoad2D::BodyForce(load))
}

fn parse_self_weight_load(parts: &[&str]) -> Result<ElementLoad2D, String> {
    let (element_id, x_acceleration, y_acceleration) = match parts {
        [_, element_id, coordinate_system, x_acceleration, y_acceleration] => {
            if !coordinate_system.eq_ignore_ascii_case("global") {
                return Err(format!("unknown self-weight coordinate system: {coordinate_system}"));
            }

            (
                parse_usize(element_id, "element ID")?,
                parse_f64(x_acceleration, "self-weight x acceleration")?,
                parse_f64(y_acceleration, "self-weight y acceleration")?,
            )
        }
        [_, element_id, x_acceleration, y_acceleration] => (
            parse_usize(element_id, "element ID")?,
            parse_f64(x_acceleration, "self-weight x acceleration")?,
            parse_f64(y_acceleration, "self-weight y acceleration")?,
        ),
        _ => return Err("expected: self_weight ELEMENT_ID global AX AY".to_owned()),
    };
    let load = SelfWeight2D::new(element_id, x_acceleration, y_acceleration).map_err(|error| error.to_string())?;

    Ok(ElementLoad2D::SelfWeight(load))
}

fn parse_load_coordinate_system(value: &str) -> Result<LoadCoordinateSystem2D, String> {
    match value.to_ascii_lowercase().as_str() {
        "global" => Ok(LoadCoordinateSystem2D::Global),
        "local" => Ok(LoadCoordinateSystem2D::Local),
        _ => Err(format!("unknown load coordinate system: {value}")),
    }
}

fn parse_section_line(line: &str) -> Result<(usize, Section2D), String> {
    let parts: Vec<&str> = line.split_whitespace().collect();

    if parts.is_empty() {
        return Err("expected a section type".to_owned());
    }

    match parts[0].to_ascii_lowercase().as_str() {
        "truss" => {
            if parts.len() != 3 {
                return Err("expected: truss ID AREA".to_owned());
            }

            let id = parse_usize(parts[1], "section ID")?;
            let cross_section_area = parse_f64(parts[2], "cross-sectional area")?;
            let section =
                TrussSection2D::new(cross_section_area).map(Section2D::Truss).map_err(|error| error.to_string())?;

            Ok((id, section))
        }
        "beam" => {
            if parts.len() != 5 {
                return Err("expected: beam ID AREA I HEIGHT".to_owned());
            }

            let id = parse_usize(parts[1], "section ID")?;
            let cross_section_area = parse_f64(parts[2], "cross-sectional area")?;
            let second_moment_of_area = parse_f64(parts[3], "second moment of area")?;
            let section_height = parse_f64(parts[4], "section height")?;
            let section =
                BeamSection2D::new_with_section_height(cross_section_area, second_moment_of_area, section_height)
                    .map(Section2D::Beam)
                    .map_err(|error| error.to_string())?;

            Ok((id, section))
        }
        "plane_stress" | "plane-stress" | "triangle" | "triangle_t3" | "t3" => {
            if parts.len() != 3 {
                return Err("expected: plane_stress ID THICKNESS".to_owned());
            }

            let id = parse_usize(parts[1], "section ID")?;
            let thickness = parse_f64(parts[2], "thickness")?;
            let section =
                PlaneStressSection2D::new(thickness).map(Section2D::PlaneStress).map_err(|error| error.to_string())?;

            Ok((id, section))
        }
        _ => Err("unknown section type; use truss, beam, or plane_stress".to_owned()),
    }
}

fn parse_element_line(line: &str) -> Result<Element2D, String> {
    let parts: Vec<&str> = line.split_whitespace().collect();

    if parts.is_empty() {
        return Err("expected an element type".to_owned());
    }

    match parts[0].to_ascii_lowercase().as_str() {
        "truss" => {
            if parts.len() != 6 {
                return Err("expected: truss ID NODE_1 NODE_2 MATERIAL_ID SECTION_ID".to_owned());
            }

            let id = parse_usize(parts[1], "element ID")?;
            let first_node_id = parse_usize(parts[2], "first node ID")?;
            let second_node_id = parse_usize(parts[3], "second node ID")?;
            let material_id = parse_usize(parts[4], "material ID")?;
            let section_id = parse_usize(parts[5], "section ID")?;

            Truss2D::new(id, [first_node_id, second_node_id], material_id, section_id)
                .map(Element2D::Truss)
                .map_err(|error| error.to_string())
        }
        "beam" => {
            if parts.len() != 6 {
                return Err("expected: beam ID NODE_1 NODE_2 MATERIAL_ID SECTION_ID".to_owned());
            }

            let id = parse_usize(parts[1], "element ID")?;
            let first_node_id = parse_usize(parts[2], "first node ID")?;
            let second_node_id = parse_usize(parts[3], "second node ID")?;
            let material_id = parse_usize(parts[4], "material ID")?;
            let section_id = parse_usize(parts[5], "section ID")?;

            Beam2D::new(id, [first_node_id, second_node_id], material_id, section_id)
                .map(Element2D::Beam)
                .map_err(|error| error.to_string())
        }
        "triangle" | "triangle_t3" | "t3" => {
            if parts.len() != 7 {
                return Err("expected: triangle ID NODE_1 NODE_2 NODE_3 MATERIAL_ID SECTION_ID".to_owned());
            }

            let id = parse_usize(parts[1], "element ID")?;
            let first_node_id = parse_usize(parts[2], "first node ID")?;
            let second_node_id = parse_usize(parts[3], "second node ID")?;
            let third_node_id = parse_usize(parts[4], "third node ID")?;
            let material_id = parse_usize(parts[5], "material ID")?;
            let section_id = parse_usize(parts[6], "section ID")?;

            TriangleT3::new(id, [first_node_id, second_node_id, third_node_id], material_id, section_id)
                .map(Element2D::TriangleT3)
                .map_err(|error| error.to_string())
        }
        _ => Err("unknown element type; use truss, beam, or triangle".to_owned()),
    }
}

fn parse_dof(value: &str) -> Result<Dof2D, String> {
    match value.to_ascii_lowercase().as_str() {
        "ux" => Ok(Dof2D::Ux),
        "uy" => Ok(Dof2D::Uy),
        "rz" => Ok(Dof2D::Rz),
        _ => Err("unknown DOF; use Ux, Uy, or Rz".to_owned()),
    }
}

fn parse_usize(value: &str, name: &str) -> Result<usize, String> {
    value.parse::<usize>().map_err(|_| format!("invalid {name}: {value}"))
}

fn parse_f64(value: &str, name: &str) -> Result<f64, String> {
    value.parse::<f64>().map_err(|_| format!("invalid {name}: {value}"))
}

#[cfg(test)]
mod tests {
    use super::{
        LoadInput2D, parse_constraint_line, parse_dof, parse_element_line, parse_load_line, parse_material_line,
        parse_node_line, parse_section_line,
    };
    use rusty_fem::elements::Element2D;
    use rusty_fem::model::{Dof2D, ElementLoad2D, LoadCoordinateSystem2D, Section2D};

    #[test]
    fn parses_node_line() {
        let node = parse_node_line("7 1.5 -2.0").expect("valid node should be parsed");

        assert_eq!(node.id(), 7);
        assert_eq!(node.coordinates(), (1.5, -2.0));
    }

    #[test]
    fn rejects_invalid_node_line() {
        assert!(parse_node_line("7 1.5").is_err());
        assert!(parse_node_line("7 abc 2.0").is_err());
    }

    #[test]
    fn parses_constraint_line() {
        let constraint = parse_constraint_line("7 Uy 0.0").expect("valid constraint should be parsed");

        assert_eq!(constraint.node_id(), 7);
        assert_eq!(constraint.dof(), Dof2D::Uy);
        assert_eq!(constraint.displacement(), 0.0);
    }

    #[test]
    fn parses_material_line() {
        let (id, material) = parse_material_line("7 200.0 0.3 7800.0").expect("valid material should be parsed");

        assert_eq!(id, 7);
        assert_eq!(material.young_modulus(), 200.0);
        assert_eq!(material.poisson_ratio(), 0.3);
        assert_eq!(material.density(), 7800.0);
    }

    #[test]
    fn rejects_invalid_material_line() {
        assert!(parse_material_line("7 200.0 0.3").is_err());
        assert!(parse_material_line("7 200.0 0.5 7800.0").is_err());
        assert!(parse_material_line("abc 200.0 0.3 7800.0").is_err());
    }

    #[test]
    fn parses_load_line() {
        let load = parse_load_line("7 Uy -12.5").expect("valid load should be parsed");
        let LoadInput2D::Nodal(load) = load else {
            panic!("expected a nodal load");
        };

        assert_eq!(load.node_id(), 7);
        assert_eq!(load.dof(), Dof2D::Uy);
        assert_eq!(load.value(), -12.5);
    }

    #[test]
    fn parses_prefixed_nodal_load_line() {
        let load = parse_load_line("nodal 7 Uy -12.5").expect("valid load should be parsed");
        let LoadInput2D::Nodal(load) = load else {
            panic!("expected a nodal load");
        };

        assert_eq!(load.node_id(), 7);
        assert_eq!(load.dof(), Dof2D::Uy);
        assert_eq!(load.value(), -12.5);
    }

    #[test]
    fn parses_beam_uniform_load_line() {
        let load = parse_load_line("beam_uniform 10 global 0.0 -5.0").expect("valid load should be parsed");
        let LoadInput2D::Element(ElementLoad2D::BeamUniformLine(load)) = load else {
            panic!("expected a uniform beam line load");
        };

        assert_eq!(load.element_id(), 10);
        assert_eq!(load.coordinate_system(), LoadCoordinateSystem2D::Global);
        assert_eq!(load.x_component(), 0.0);
        assert_eq!(load.y_component(), -5.0);
    }

    #[test]
    fn parses_edge_traction_load_line() {
        let load = parse_load_line("edge_traction 10 2 3 local 4.0 -8.0").expect("valid load should be parsed");
        let LoadInput2D::Element(ElementLoad2D::EdgeTraction(load)) = load else {
            panic!("expected an edge traction load");
        };

        assert_eq!(load.element_id(), 10);
        assert_eq!(load.edge_node_ids(), [2, 3]);
        assert_eq!(load.coordinate_system(), LoadCoordinateSystem2D::Local);
        assert_eq!(load.x_component(), 4.0);
        assert_eq!(load.y_component(), -8.0);
    }

    #[test]
    fn parses_body_force_load_line() {
        let load = parse_load_line("body_force 10 global 4.0 -8.0").expect("valid load should be parsed");
        let LoadInput2D::Element(ElementLoad2D::BodyForce(load)) = load else {
            panic!("expected a body force load");
        };

        assert_eq!(load.element_id(), 10);
        assert_eq!(load.x_component(), 4.0);
        assert_eq!(load.y_component(), -8.0);
    }

    #[test]
    fn parses_body_force_load_line_without_global_keyword() {
        let load = parse_load_line("body_force 10 4.0 -8.0").expect("valid load should be parsed");
        let LoadInput2D::Element(ElementLoad2D::BodyForce(load)) = load else {
            panic!("expected a body force load");
        };

        assert_eq!(load.element_id(), 10);
        assert_eq!(load.x_component(), 4.0);
        assert_eq!(load.y_component(), -8.0);
    }

    #[test]
    fn parses_self_weight_load_line() {
        let load = parse_load_line("self_weight 10 global 0.0 -9.81").expect("valid load should be parsed");
        let LoadInput2D::Element(ElementLoad2D::SelfWeight(load)) = load else {
            panic!("expected a self-weight load");
        };

        assert_eq!(load.element_id(), 10);
        assert_eq!(load.x_acceleration(), 0.0);
        assert_eq!(load.y_acceleration(), -9.81);
    }

    #[test]
    fn parses_self_weight_load_line_without_global_keyword() {
        let load = parse_load_line("self_weight 10 0.0 -9.81").expect("valid load should be parsed");
        let LoadInput2D::Element(ElementLoad2D::SelfWeight(load)) = load else {
            panic!("expected a self-weight load");
        };

        assert_eq!(load.element_id(), 10);
        assert_eq!(load.x_acceleration(), 0.0);
        assert_eq!(load.y_acceleration(), -9.81);
    }

    #[test]
    fn rejects_invalid_load_line() {
        assert!(parse_load_line("7 Uy").is_err());
        assert!(parse_load_line("7 Rz NaN").is_err());
        assert!(parse_load_line("7 Unknown 10.0").is_err());
        assert!(parse_load_line("beam_uniform 10 global 0").is_err());
        assert!(parse_load_line("beam_uniform 10 rotated 0 -5").is_err());
        assert!(parse_load_line("edge_traction 10 2 3 global 0").is_err());
        assert!(parse_load_line("edge_traction 10 2 2 global 0 -5").is_err());
        assert!(parse_load_line("edge_traction 10 2 3 rotated 0 -5").is_err());
        assert!(parse_load_line("body_force 10 global 0").is_err());
        assert!(parse_load_line("body_force 10 local 0 -5").is_err());
        assert!(parse_load_line("self_weight 10 global 0").is_err());
        assert!(parse_load_line("self_weight 10 local 0 -9.81").is_err());
    }

    #[test]
    fn parses_dofs_case_insensitively() {
        assert_eq!(parse_dof("ux"), Ok(Dof2D::Ux));
        assert_eq!(parse_dof("Uy"), Ok(Dof2D::Uy));
        assert_eq!(parse_dof("RZ"), Ok(Dof2D::Rz));
    }

    #[test]
    fn parses_section_line() {
        let truss = parse_section_line("truss 10 0.01").expect("valid truss section should be parsed");
        let beam = parse_section_line("beam 15 0.02 0.001 0.1").expect("valid beam section should be parsed");
        let plane_stress =
            parse_section_line("plane_stress 20 0.1").expect("valid plane-stress section should be parsed");

        assert!(matches!(truss, (10, Section2D::Truss(_))));
        assert!(matches!(
            beam,
            (15, Section2D::Beam(section)) if section.section_height() == Some(0.1)
        ));
        assert!(matches!(plane_stress, (20, Section2D::PlaneStress(_))));
    }

    #[test]
    fn rejects_invalid_section_line() {
        assert!(parse_section_line("beam 10 0.02").is_err());
        assert!(parse_section_line("beam 10 0.02 0.001").is_err());
        assert!(parse_section_line("truss 10").is_err());
        assert!(parse_section_line("plane_stress 10").is_err());
        assert!(parse_section_line("hexagon 10 0.1").is_err());
    }

    #[test]
    fn parses_element_line() {
        let truss = parse_element_line("truss 10 1 2 7 100").expect("valid truss should be parsed");
        let beam = parse_element_line("beam 15 1 2 8 200").expect("valid beam should be parsed");
        let triangle = parse_element_line("triangle 20 1 2 3 9 300").expect("valid triangle should be parsed");

        assert!(
            matches!(truss, Element2D::Truss(element) if element.material_id() == 7 && element.section_id() == 100)
        );
        assert!(matches!(beam, Element2D::Beam(element) if element.material_id() == 8 && element.section_id() == 200));
        assert!(
            matches!(triangle, Element2D::TriangleT3(element) if element.material_id() == 9 && element.section_id() == 300)
        );
    }

    #[test]
    fn rejects_invalid_element_line() {
        assert!(parse_element_line("beam 10 1").is_err());
        assert!(parse_element_line("beam 10 1 2").is_err());
        assert!(parse_element_line("beam 10 1 2 7").is_err());
        assert!(parse_element_line("truss 10 1 2").is_err());
        assert!(parse_element_line("triangle 10 1 2 3 7").is_err());
        assert!(parse_element_line("hexagon 10 1 2").is_err());
    }
}
