// Line shader: simple pass-through for wireframe overlay rendering.

struct Uniforms {
    view_proj: mat4x4<f32>,
    light_dir: vec4<f32>,
};

@group(0) @binding(0) var<uniform> uniforms: Uniforms;

struct LineInput {
    @location(0) position: vec3<f32>,
    @location(1) color: vec4<f32>,
};

struct LineOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) color: vec4<f32>,
};

@vertex
fn vs_main(in: LineInput) -> LineOutput {
    var out: LineOutput;
    out.clip_position = uniforms.view_proj * vec4<f32>(in.position, 1.0);
    out.color = in.color;
    return out;
}

@fragment
fn fs_main(in: LineOutput) -> @location(0) vec4<f32> {
    return in.color;
}
