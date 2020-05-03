
void main_1() {
    vec2 rg = iMouse.xy / iResolution.xy;
    gl_FragColor = vec4(rg, 0.2, 1);
}

void main_2() {
    vec2 uv = gl_FragCoord.xy / iResolution.xy;

    float y = pow(uv.x, 5.0);
    
}

void main() {
    main_1();
}
