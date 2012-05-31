attribute vec4 position;
attribute vec2 textureCoordinates;

varying vec2 v_textureCoordinates;

void main()
{
    v_textureCoordinates = textureCoordinates;
    gl_Position = agi_viewportOrthographic * position;
}